#!/bin/bash
# =============================================================================
# make_ws_and_datacard.sh
#
# Full Run3 VBF H→inv pipeline with automatic environment switching:
#
#   Step 1  NLO ratio histograms  (derive_run3_nlo_ratios.py)   — pyRAT/LCG
#   Step 2  Signal + MC bkg WS    (makeSignalAndMCBackgroundWSRun3.C) — CMSSW
#   Step 3  Parametric W/Z TF WS  (makeWS_percategoryRun3.C)    — CMSSW
#   Step 4  Datacard              (make_run3_datacard.py)        — pyRAT/LCG
#   Step 5  Combine cards         (combineCards.py)             — CMSSW (multi-cat only)
#   Step 6  Expected limit        (combine AsymptoticLimits)    — CMSSW
#
# LCG_104 (loaded by setup_environment.sh) and CMSSW are incompatible in the
# same shell: they conflict on ROOTSYS, LD_LIBRARY_PATH, and PYTHONPATH.
# Each step therefore runs in a fresh subprocess (via env -i + bash) with only
# its required environment sourced.  No CMSSW pre-activation is needed. 
#
# Usage (from /vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT):
#   bash makeHinvWS/make_ws_and_datacard.sh
#
# Optional env vars:
#   CATS="VTR MTR"                                         (default: VTR)
#   VAR="Mjj"                                              (default: Mjj;  also: SignalScore)
#   YEAR="Run3Summer22_to_Run3Summer23BPix"                (default: Run3Summer22_to_Run3Summer23BPix)
#   CAMPAIGNS="Run3Summer22 Run3Summer22EE Run3Summer23"   (default: all four Run3 campaigns)
#   SKIP_LIMITS=1     skip Step 6 (just build workspaces + datacard)
#   OBSERVED=1        run combine on observed data instead of blinded
# =============================================================================

set -euo pipefail

# Parse key=value positional arguments (supports both `bash script.sh VAR=X`
# and `source script.sh VAR=X` invocation styles)
for _arg in "$@"; do
    case "${_arg}" in
        CATS=*)      CATS="${_arg#CATS=}" ;;
        VAR=*)       VAR="${_arg#VAR=}" ;;
        YEAR=*)      YEAR="${_arg#YEAR=}" ;;
        CAMPAIGNS=*) CAMPAIGNS="${_arg#CAMPAIGNS=}" ;;
        SKIP_LIMITS=*) SKIP_LIMITS="${_arg#SKIP_LIMITS=}" ;;
        OBSERVED=*)    OBSERVED="${_arg#OBSERVED=}" ;;
    esac
done
unset _arg

# ---- configuration ----------------------------------------------------------
CATS="${CATS:-VTR}"
VAR="${VAR:-Mjj}"
[ "${VAR}" = "SignalScore" ] && CLASSIFIER="true" || CLASSIFIER="false"
YEAR="${YEAR:-Run3Summer22_to_Run3Summer23BPix}"
CAMPAIGNS="${CAMPAIGNS:-Run3Summer22 Run3Summer22EE Run3Summer23 Run3Summer23BPix}"

PYRAT="/vols/cms/tt1020/HiggsInvisible/pyRAT/pyRAT"
MKHWS="${PYRAT}/makeHinvWS"
DATACARDS="${PYRAT}/datacards"
CMSSW_PATH="/vols/cms/tt1020/Combine/CMSSW_14_1_0_pre4"

ROOT_CMD="root -l -b -q"
RECOMPILE="++"   # ++ = always recompile; + = skip if .so is current

# Capture HOME now — env -i subshells need it passed explicitly
_HOME="${HOME}"

# =============================================================================
# Environment helpers
#
# run_in_pyrat  "cmd"  — strips env, sources LCG_104 + kraken venv, runs cmd
# run_in_cmssw  "cmd"  — strips env, sources CMSSW_14_1_0_pre4, runs cmd
#
# Using env -i means neither environment leaks into the other: the three
# variables that cause conflicts (ROOTSYS, LD_LIBRARY_PATH, PYTHONPATH) are
# absent from the subprocess before the correct env is sourced.
# =============================================================================

run_in_pyrat() {
    # $1 = command string to run inside pyRAT/LCG environment
    env -i HOME="${_HOME}" bash -c "
        export PATH=/usr/local/bin:/usr/bin:/bin
        source '${PYRAT}/setup_environment.sh'
        $1
    "
}

run_in_cmssw() {
    # $1 = command string to run inside CMSSW environment
    env -i HOME="${_HOME}" bash -c "
        export PATH=/usr/local/bin:/usr/bin:/bin
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        cd '${CMSSW_PATH}/src' && eval \$(scramv1 runtime -sh) && cd '${PYRAT}'
        $1
    "
}

# =============================================================================

echo ""
echo "======================================================================"
echo "  Run3 VBF H→inv pipeline"
echo "  Categories : ${CATS}"
echo "  Variable   : ${VAR}"
echo "  Year       : ${YEAR}"
echo "======================================================================"

# ---- per-category steps 1–4 -------------------------------------------------
for CAT in ${CATS}; do

    echo ""
    echo "----------------------------------------------------------------------"
    echo "  Category: ${CAT}"
    echo "----------------------------------------------------------------------"

    LABEL="${CAT}_${YEAR}"

    # ------------------------------------------------------------------
    # Step 1: NLO ratio histograms  [pyRAT/LCG]
    # Reads SR VBF_shapes.root → run3_nlo_sf_*.root (needed by Step 3)
    # ------------------------------------------------------------------
    echo ""
    echo "[Step 1/${CAT}]  NLO ratio histograms  (pyRAT env)"
    run_in_pyrat "python3 '${MKHWS}/derive_run3_nlo_ratios.py' \
        --year  '${YEAR}' \
        --cat   '${CAT}'  \
        --var   '${VAR}'  \
        --output-dir '${MKHWS}'"

    # ------------------------------------------------------------------
    # Step 2: Signal + MC background workspace  [CMSSW]
    # makeSignalAndMCBackgroundWSRun3.C → signal_mc_bkgs_ws_*.root
    # ------------------------------------------------------------------
    echo ""
    echo "[Step 2/${CAT}]  Signal + MC background workspace  (CMSSW env)"
    run_in_cmssw "cd '${MKHWS}' && ${ROOT_CMD} \
        'makeSignalAndMCBackgroundWSRun3.C${RECOMPILE}(\"${YEAR}\",\"${CAT}\",${CLASSIFIER},true,true)'"

    # ------------------------------------------------------------------
    # Step 3: Parametric W/Z TF workspace  [CMSSW]
    # makeWS_percategoryRun3.C → param_ws_*.root
    # MUST come after Step 1 (reads run3_nlo_sf_*.root from makeHinvWS/).
    # ------------------------------------------------------------------
    echo ""
    echo "[Step 3/${CAT}]  Parametric W/Z TF workspace  (CMSSW env)"
    run_in_cmssw "cd '${MKHWS}' && ${ROOT_CMD} \
        'makeWS_percategoryRun3.C${RECOMPILE}(\"${YEAR}\",\"${CAT}\",${CLASSIFIER},true,true)'"

    # ------------------------------------------------------------------
    # Step 4: Datacard  [pyRAT/LCG]
    # Copies workspaces + DD files to datacards/{LABEL}/ and writes .txt
    # ------------------------------------------------------------------
    echo ""
    echo "[Step 4/${CAT}]  Datacard  (pyRAT env)"
    run_in_pyrat "python3 '${MKHWS}/make_run3_datacard.py' \
        --cat       '${CAT}'     \
        --campaigns ${CAMPAIGNS} \
        --var       '${VAR}'     \
        --label     '${LABEL}'"

done

# ---- Step 5: combine cards (multi-category, CMSSW) --------------------------
NCATS=$(echo "${CATS}" | wc -w)

if [ "${NCATS}" -gt 1 ]; then
    echo ""
    echo "[Step 5]  Combining per-category datacards  (CMSSW env)"
    CARD_ARGS=""
    for CAT in ${CATS}; do
        LABEL="${CAT}_${YEAR}"
        CARD_ARGS="${CARD_ARGS} ${CAT}_Run3=${DATACARDS}/${LABEL}/card_${LABEL}.txt"
    done
    COMBINED="${DATACARDS}/run3_combined_${YEAR}.txt"
    run_in_cmssw "combineCards.py ${CARD_ARGS} > '${COMBINED}'"
    echo "  Combined card: ${COMBINED}"
    FIT_CARD="${COMBINED}"
else
    CAT="${CATS}"
    LABEL="${CAT}_${YEAR}"
    FIT_CARD="${DATACARDS}/${LABEL}/card_${LABEL}.txt"
fi

# ---- Step 6: expected limit (CMSSW) -----------------------------------------
if [ "${SKIP_LIMITS:-0}" = "1" ]; then
    echo ""
    echo "[Step 6]  Skipped (SKIP_LIMITS=1)"
else
    BLIND_FLAG="--run blind -t -1"
    [ "${OBSERVED:-0}" = "1" ] && BLIND_FLAG=""

    echo ""
    echo "[Step 6]  AsymptoticLimits  (CMSSW env)"
    LOG="${DATACARDS}/combine_run3_${YEAR}.log"
    run_in_cmssw "cd '${DATACARDS}' && \
        combine -M AsymptoticLimits \
            -d '${FIT_CARD}' \
            -m 125 \
            ${BLIND_FLAG} \
            -n 'Run3_${YEAR}' \
            2>&1 | tee '${LOG}'"

    echo ""
    echo "  Results:"
    grep -E "Expected|Observed" "${LOG}" || true
    echo ""
    echo "  Output: higgsCombineRun3_${YEAR}.AsymptoticLimits.mH125.root"
fi

echo ""
echo "======================================================================"
echo "  Pipeline complete."
echo "======================================================================"
