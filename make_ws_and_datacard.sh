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
#   Step 5  Combine cards         (combineCards.py)             — CMSSW (multi-cat or multi-year)
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
#   YEAR=""                (default: Run3Summer22_to_Run3Summer23BPix)
#   YEARS="Run3Summer22_to_Run3Summer23BPix RunIII2024Summer24"
#                          (space-separated list; overrides YEAR when set)
#   CAMPAIGNS="Run3Summer22 Run3Summer22EE Run3Summer23 Run3Summer23BPix"
#                          (explicit campaign list; only honoured for single-year runs)
#   SKIP_LIMITS=1     skip Step 6 (just build workspaces + datacard)
#   OBSERVED=1        run combine on observed data instead of blinded
#   JEC_MODE="none"   JEC treatment: "none" (off), "dummy" (flat Run2-derived mean, default), "maxdummy" (flat Run2-derived max), "real"
#   PHOTON=1          also build the γ+jets CR workspace (makePhotonWSRun3.C) and include it in the datacard
#   NO_HF_NOISE=1     exclude the HF noise data-driven estimate from the SR and datacard
# =============================================================================

set -euo pipefail

# Parse key=value positional arguments (supports both `bash script.sh VAR=X`
# and `source script.sh VAR=X` invocation styles)
_CAMPAIGNS_EXPLICIT=0
for _arg in "$@"; do
    case "${_arg}" in
        CATS=*)        CATS="${_arg#CATS=}" ;;
        VAR=*)         VAR="${_arg#VAR=}" ;;
        YEAR=*)        YEAR="${_arg#YEAR=}" ;;
        YEARS=*)       YEARS="${_arg#YEARS=}" ;;
        CAMPAIGNS=*)   CAMPAIGNS="${_arg#CAMPAIGNS=}"; _CAMPAIGNS_EXPLICIT=1 ;;
        SKIP_LIMITS=*)  SKIP_LIMITS="${_arg#SKIP_LIMITS=}" ;;
        OBSERVED=*)     OBSERVED="${_arg#OBSERVED=}" ;;
        JEC_MODE=*)     JEC_MODE="${_arg#JEC_MODE=}" ;;
        PHOTON=*)       PHOTON="${_arg#PHOTON=}" ;;
        NO_HF_NOISE=*)  NO_HF_NOISE="${_arg#NO_HF_NOISE=}" ;;
    esac
done
unset _arg

# ---- configuration ----------------------------------------------------------
CATS="${CATS:-VTR}"
VAR="${VAR:-Mjj}"
[ "${VAR}" = "SignalScore" ] && CLASSIFIER="true" || CLASSIFIER="false"
YEAR="${YEAR:-Run3Summer22_to_Run3Summer23BPix}"
YEARS="${YEARS:-${YEAR}}"
NYEARS=$(echo "${YEARS}" | wc -w)
NCATS=$(echo "${CATS}" | wc -w)

# campaigns_for_year <year>
# Auto-derives the CAMPAIGNS list for a given YEAR string.
# For single-year runs an explicit CAMPAIGNS= override is honoured.
campaigns_for_year() {
    local _year="$1"
    if [ "${_CAMPAIGNS_EXPLICIT}" = "1" ] && [ "${NYEARS}" -eq 1 ]; then
        echo "${CAMPAIGNS}"
        return
    fi
    case "${_year}" in
        Run3Summer22_to_Run3Summer23BPix)
            echo "Run3Summer22 Run3Summer22EE Run3Summer23 Run3Summer23BPix" ;;
        *)
            echo "${_year}" ;;
    esac
}

JEC_MODE="${JEC_MODE:-dummy}"   # "none" / "dummy" / "maxdummy" / "real"
PHOTON="${PHOTON:-0}"           # "0" off  "1" build photon CR workspace + include in datacard
NO_HF_NOISE="${NO_HF_NOISE:-0}" # "0" include HF noise DD estimate  "1" exclude it
# Translate maxdummy -> dummy for the C++ scripts; the difference is only in the input ROOT files
_CPP_JEC_MODE="${JEC_MODE}"
[ "${JEC_MODE}" = "maxdummy" ] && _CPP_JEC_MODE="dummy"

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
echo "  Years      : ${YEARS}"
echo "  Photon CR  : ${PHOTON}"
echo "======================================================================"

# ---- Step 1b: regenerate JEC dummy uncertainty files ------------------------
# Runs once (not per-category or per-year) when JEC_MODE is dummy or maxdummy.
# Always regenerates for ALL known Run3 years so no histogram is lost.
if [ "${JEC_MODE}" = "dummy" ] || [ "${JEC_MODE}" = "maxdummy" ]; then
    _MAX_FLAG=""
    [ "${JEC_MODE}" = "maxdummy" ] && _MAX_FLAG="--use_max"
    echo ""
    echo "[Step 1b]  JEC dummy uncertainty files  (pyRAT env)  mode=${JEC_MODE}"
    run_in_pyrat "python3 '${MKHWS}/make_run3_jec_uncs.py' \
        --run2_shape '${MKHWS}/vbf_shape_jes_uncs.root' \
        --run2_tf    '${MKHWS}/vbf_jes_jer_tf_uncs.root' \
        --out_shape  '${MKHWS}/vbf_shape_jes_uncs_run3.root' \
        --out_tf     '${MKHWS}/vbf_jes_jer_tf_uncs_run3.root' \
        ${_MAX_FLAG}"
fi

# ---- per-year, per-category steps 1–4 --------------------------------------
for _YEAR in ${YEARS}; do
    _CAMPAIGNS="$(campaigns_for_year "${_YEAR}")"

    echo ""
    echo "======================================================================"
    echo "  Year: ${_YEAR}  |  Campaigns: ${_CAMPAIGNS}"
    echo "======================================================================"

    for CAT in ${CATS}; do

        echo ""
        echo "----------------------------------------------------------------------"
        echo "  Category: ${CAT}"
        echo "----------------------------------------------------------------------"

        LABEL="${CAT}_${_YEAR}"

        # ------------------------------------------------------------------
        # Step 1: NLO ratio histograms  [pyRAT/LCG]
        # Reads SR VBF_shapes.root → run3_nlo_sf_*.root (needed by Step 3)
        # ------------------------------------------------------------------
        echo ""
        echo "[Step 1/${CAT}/${_YEAR}]  NLO ratio histograms  (pyRAT env)"
        run_in_pyrat "python3 '${MKHWS}/derive_run3_nlo_ratios.py' \
            --year  '${_YEAR}' \
            --cat   '${CAT}'   \
            --var   '${VAR}'   \
            --output-dir '${MKHWS}'"

        # ------------------------------------------------------------------
        # Step 2: Signal + MC background workspace  [CMSSW]
        # makeSignalAndMCBackgroundWSRun3.C → signal_mc_bkgs_ws_*.root
        # ------------------------------------------------------------------
        echo ""
        echo "[Step 2/${CAT}/${_YEAR}]  Signal + MC background workspace  (CMSSW env)"
        run_in_cmssw "cd '${MKHWS}' && ${ROOT_CMD} \
            'makeSignalAndMCBackgroundWSRun3.C${RECOMPILE}(\"${_YEAR}\",\"${CAT}\",${CLASSIFIER},\"${_CPP_JEC_MODE}\",true)'"

        # ------------------------------------------------------------------
        # Step 2b: Photon CR workspace  [CMSSW]  (only when PHOTON=1)
        # makePhotonWSRun3.C reads singlePhoton_CR shapes and SR QCD_Znunu shapes,
        # then builds parametric γ+jets PDFs linked to the shared Z SR free parameters.
        # ------------------------------------------------------------------
        if [ "${PHOTON}" = "1" ]; then
            echo ""
            echo "[Step 2b/${CAT}/${_YEAR}]  Photon CR workspace  (CMSSW env)"
            run_in_cmssw "cd '${MKHWS}' && ${ROOT_CMD} \
                'makePhotonWSRun3.C${RECOMPILE}(\"${_YEAR}\",\"${CAT}\",${CLASSIFIER},\"${_CPP_JEC_MODE}\",true)'"
        fi

        # ------------------------------------------------------------------
        # Step 3: Parametric W/Z TF workspace  [CMSSW]
        # makeWS_percategoryRun3.C → param_ws_*.root
        # MUST come after Step 1 (reads run3_nlo_sf_*.root from makeHinvWS/).
        # ------------------------------------------------------------------
        echo ""
        echo "[Step 3/${CAT}/${_YEAR}]  Parametric W/Z TF workspace  (CMSSW env)"
        run_in_cmssw "cd '${MKHWS}' && ${ROOT_CMD} \
            'makeWS_percategoryRun3.C${RECOMPILE}(\"${_YEAR}\",\"${CAT}\",${CLASSIFIER},\"${_CPP_JEC_MODE}\",true)'"

        # ------------------------------------------------------------------
        # Step 4: Datacard  [pyRAT/LCG]
        # Copies workspaces + DD files to datacards/{LABEL}/ and writes .txt
        # ------------------------------------------------------------------
        echo ""
        echo "[Step 4/${CAT}/${_YEAR}]  Datacard  (pyRAT env)"
        _PHOTON_FLAG=""
        [ "${PHOTON}" = "1" ] && _PHOTON_FLAG="--photon"
        _NO_HF_NOISE_FLAG=""
        [ "${NO_HF_NOISE}" = "1" ] && _NO_HF_NOISE_FLAG="--no-hf-noise"
        run_in_pyrat "python3 '${MKHWS}/make_run3_datacard.py' \
            --cat       '${CAT}'       \
            --campaigns ${_CAMPAIGNS}  \
            --var       '${VAR}'       \
            --label     '${LABEL}'     \
            ${_PHOTON_FLAG} ${_NO_HF_NOISE_FLAG}"

    done
done

# ---- Step 5: combine cards (any time more than one card was built) ----------
NTOTAL=$((NYEARS * NCATS))

if [ "${NTOTAL}" -gt 1 ]; then
    echo ""
    echo "[Step 5]  Combining per-category/per-year datacards  (CMSSW env)"
    CARD_ARGS=""
    for _YEAR in ${YEARS}; do
        for CAT in ${CATS}; do
            LABEL="${CAT}_${_YEAR}"
            CARD_ARGS="${CARD_ARGS} ${CAT}_${_YEAR}=${DATACARDS}/${LABEL}/card_${LABEL}.txt"
        done
    done
    _YEARS_SLUG=$(echo "${YEARS}" | tr ' ' '_plus_')
    COMBINED="${DATACARDS}/run3_combined_${_YEARS_SLUG}.txt"
    run_in_cmssw "combineCards.py ${CARD_ARGS} > '${COMBINED}'"
    echo "  Combined card: ${COMBINED}"
    FIT_CARD="${COMBINED}"
else
    _YEAR="${YEARS}"
    CAT="${CATS}"
    LABEL="${CAT}_${_YEAR}"
    FIT_CARD="${DATACARDS}/${LABEL}/card_${LABEL}.txt"
fi

# ---- Step 6: expected limit (CMSSW) -----------------------------------------
if [ "${SKIP_LIMITS:-0}" = "1" ]; then
    echo ""
    echo "[Step 6]  Skipped (SKIP_LIMITS=1)"
else
    BLIND_FLAG="--run blind -t -1"
    [ "${OBSERVED:-0}" = "1" ] && BLIND_FLAG=""

    _YEARS_SLUG=$(echo "${YEARS}" | tr ' ' '_plus_')
    echo ""
    echo "[Step 6]  AsymptoticLimits  (CMSSW env)"
    LOG="${DATACARDS}/combine_run3_${_YEARS_SLUG}.log"
    run_in_cmssw "cd '${DATACARDS}' && \
        combine -M AsymptoticLimits \
            -d '${FIT_CARD}' \
            -m 125 \
            ${BLIND_FLAG} \
            -n 'Run3_${_YEARS_SLUG}' \
            2>&1 | tee '${LOG}'"

    echo ""
    echo "  Results:"
    grep -E "Expected|Observed" "${LOG}" || true
    echo ""
    echo "  Output: higgsCombineRun3_${_YEARS_SLUG}.AsymptoticLimits.mH125.root"
fi

echo ""
echo "======================================================================"
echo "  Pipeline complete."
echo "======================================================================"
