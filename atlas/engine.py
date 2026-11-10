"""
Growth-kinetics engine used to compute atlas phenotypes.

Copied unchanged from the Bioprocess Analyzer
(github.com/anupama29k/bioprocessanalyzer, modules/calcs.py @ 2439bf9) so the atlas
builds without installing the analyzer app. Keep in sync when the engine changes:
an engine change is a new atlas release.
"""
import numpy as np
from scipy.stats import linregress

ENGINE_SOURCE = "bioprocessanalyzer@2439bf9"

def _detect_growth_changepoints(t, ln_od, min_seg=None, penalty_factor=None,
                                 _depth=0, _max_depth=3):
    """Detect changepoints in ln(OD) where the growth rate shifts.

    Uses penalised piecewise-linear regression: tries every candidate split
    and keeps those that reduce total RSS by more than a noise-based penalty.
    Returns a sorted list of changepoint indices (excluding endpoints).

    Guards against over-segmentation via:
    - Adaptive min_seg based on data density (at least 8 points or 10 %)
    - Adaptive penalty: tighter for short/noisy data, relaxed for long runs
    - Minimum slope ratio: slopes must differ by >= 40 % relative
    - Max recursion depth of 3 (max 4 segments)
    """
    n = len(t)
    if min_seg is None:
        min_seg = max(8, n // 10)
    if penalty_factor is None:
        # Longer datasets are more likely to contain true phase shifts,
        # so use a lower penalty. Short/noisy data gets a higher penalty
        # to prevent noise-driven splits.
        penalty_factor = 0.20 if n >= 40 else 0.40
    if n < 2 * min_seg or _depth >= _max_depth:
        return []

    # Single-segment baseline RSS
    sl0, ic0, _, _, _ = linregress(t, ln_od)
    rss0 = float(np.sum((ln_od - (sl0 * t + ic0)) ** 2))
    noise_var = rss0 / max(n - 2, 1)
    penalty = penalty_factor * n * noise_var

    best_gain, best_cp = -np.inf, None

    for cp in range(min_seg, n - min_seg):
        sl_l, ic_l, _, _, _ = linregress(t[:cp], ln_od[:cp])
        rss_l = float(np.sum((ln_od[:cp] - (sl_l * t[:cp] + ic_l)) ** 2))
        sl_r, ic_r, _, _, _ = linregress(t[cp:], ln_od[cp:])
        rss_r = float(np.sum((ln_od[cp:] - (sl_r * t[cp:] + ic_r)) ** 2))

        gain = rss0 - (rss_l + rss_r) - penalty
        # Slopes must differ by >= 40 % relative to the steeper one
        max_slope = max(abs(sl_l), abs(sl_r), 1e-9)
        slope_ratio = abs(sl_l - sl_r) / max_slope
        if gain > best_gain and slope_ratio > 0.40:
            best_gain = gain
            best_cp = cp

    if best_cp is not None and best_gain > 0:
        left_cps  = _detect_growth_changepoints(
            t[:best_cp], ln_od[:best_cp], min_seg, penalty_factor,
            _depth + 1, _max_depth)
        right_cps = _detect_growth_changepoints(
            t[best_cp:], ln_od[best_cp:], min_seg, penalty_factor,
            _depth + 1, _max_depth)
        right_cps = [cp + best_cp for cp in right_cps]
        return left_cps + [best_cp] + right_cps

    return []


MAX_LN_DROP = np.log(0.75)


def _has_plateau(t, ln_od):
    """True if the last readings grow much slower than the fastest stretch."""
    if len(t) < 3:
        return False
    rates = np.diff(ln_od) / np.maximum(np.diff(t), 1e-9)
    peak = np.max(rates)
    if peak <= 0:
        return True
    return bool(rates[-1] < 0.25 * peak)


def _fit_phase(t, od, ln_od, phase_start, phase_end, min_pts=5, plateaued=True):
    """Fit a single growth phase using sliding-window regression.

    Returns dict with mu, r2, etc. for this phase, or None if no fit.
    """
    n_ph = phase_end - phase_start
    t_ph = t[phase_start:phase_end]
    od_ph = od[phase_start:phase_end]
    ln_ph = ln_od[phase_start:phase_end]
    max_od_ph = float(np.max(od_ph))

    if n_ph < min_pts:
        return None

    span_ph = float(t_ph[-1] - t_ph[0])
    # 15 % of the phase, capped at 6 h so long stationary tails (fed-batch,
    # extended runs) do not demand a window longer than the exponential phase.
    min_span = max(2.0, min(span_ph * 0.15, 6.0))
    best_score, best = -np.inf, None

    for i in range(n_ph - min_pts + 1):
        if i > 0 and od_ph[i] < od_ph[0] * 0.40:
            continue
        for j in range(i + min_pts, n_ph + 1):
            span = t_ph[j-1] - t_ph[i]
            if span < min_span:
                continue
            if plateaued and od_ph[i:j].max() > max_od_ph * 0.92:
                continue
            # A true exponential window cannot contain a large OD drop.
            # Drops >25 % between consecutive readings are near-detection-
            # limit noise or sampling artefacts; including them inflates μ.
            if np.min(np.diff(ln_ph[i:j])) < MAX_LN_DROP:
                continue
            sl, ic, r, _, _ = linregress(t_ph[i:j], ln_ph[i:j])
            if sl <= 0:
                continue
            r2 = r ** 2
            # Noise-penalised scoring: R²² down-weights windows with
            # marginal fit; rstd penalty discourages absorbing noisy
            # low-OD tails that bias the slope downward.
            resids = ln_ph[i:j] - (sl * t_ph[i:j] + ic)
            rstd = float(np.std(resids))
            score = sl * r2 ** 2 * np.sqrt(span) / (1.0 + 4.0 * rstd)
            if score > best_score:
                best_score = score
                best = (i, j, sl, ic, r2)

    if best is None:
        return None

    i0, j, sl, ic, r2 = best

    # Post-selection trim
    _resids = ln_ph[i0:j] - (sl * t_ph[i0:j] + ic)
    _rstd = max(float(np.std(_resids)), 0.10)
    if _resids[0] < -1.5 * _rstd and (j - i0 - 1) >= min_pts:
        sl2, ic2, r2_2, _, _ = linregress(t_ph[i0+1:j], ln_ph[i0+1:j])
        if sl2 > 0 and r2_2 >= r2 * 0.97:
            i0, sl, ic, r2 = i0 + 1, float(sl2), float(ic2), float(r2_2)

    _, _, _, _, mu_se = linregress(t_ph[i0:j], ln_ph[i0:j])
    n_exp = j - i0

    return {
        "mu":       float(sl),
        "mu_se":    float(mu_se),
        "r2":       float(r2),
        "n_exp":    int(n_exp),
        "td_min":   float(np.log(2) / sl * 60),
        "t_start":  float(t_ph[0]),
        "t_end":    float(t_ph[-1]),
        "exp_t":    t_ph[i0:j],
        "bgr_sl":   float(sl),
        "bgr_ic":   float(ic),
        "i0_global": phase_start + i0,
        "j_global":  phase_start + j,
    }


def compute_growth_kinetics(time_h, od, dcw_factor=0.4):
    t  = np.asarray(time_h, float)
    od = np.asarray(od, float)
    dcw    = od * dcw_factor
    max_od = float(np.max(od))
    res    = {"dcw": dcw, "max_od": max_od, "max_dcw": float(np.max(dcw))}

    n          = len(t)
    total_span = float(t[-1] - t[0])
    ln_od      = np.log(np.clip(od, 1e-9, None))
    # Bioreactor runs are often sparsely sampled (5–8 OD points); allow a
    # 3-point exponential window there, flagged by a lower quality score.
    MIN_PTS    = 5 if n >= 10 else max(3, n - 2)
    # Exclude near-maximum readings only if the whole culture plateaued; runs
    # harvested mid-exponential have their highest OD still on the line.
    plateaued  = _has_plateau(t, ln_od)

    # ── Diauxic / multi-phase detection ──────────────────────────────────────
    # Detect changepoints where growth rate shifts significantly, then fit
    # each segment independently.  The primary μ is the fastest phase.
    changepoints = _detect_growth_changepoints(t, ln_od)
    boundaries = [0] + changepoints + [n]

    phases = []
    for seg_idx in range(len(boundaries) - 1):
        seg_start = boundaries[seg_idx]
        seg_end   = boundaries[seg_idx + 1]
        phase = _fit_phase(t, od, ln_od, seg_start, seg_end, MIN_PTS, plateaued)
        if phase is not None:
            phase["phase_idx"] = seg_idx
            phases.append(phase)

    # If segmented fitting yielded no phases, fall back to full-dataset fit
    if not phases:
        full_phase = _fit_phase(t, od, ln_od, 0, n, MIN_PTS, plateaued)
        if full_phase is not None:
            full_phase["phase_idx"] = 0
            phases.append(full_phase)
        changepoints = []

    # Always also fit the full dataset: changepoint boundaries can clip
    # valid exponential data, so the unconstrained fit is a useful
    # competitor.  It's added as a hidden candidate — not a reported phase.
    full_fit = _fit_phase(t, od, ln_od, 0, n, MIN_PTS, plateaued) if changepoints else None

    res["phases"] = phases
    res["changepoints"] = [float(t[cp]) for cp in changepoints]
    is_diauxic = len(phases) >= 2
    res["is_diauxic"] = is_diauxic

    # ── Select primary phase (fastest μ) for backward-compatible keys ────────
    if phases:
        primary = max(phases, key=lambda p: p["mu"])
        # If the full-dataset fit found a higher-R² result at similar μ,
        # prefer it (changepoint boundary may have clipped good data).
        if full_fit is not None and full_fit["r2"] > primary["r2"] + 0.005:
            # Only override if the mu values are reasonably close (within 30%)
            mu_ratio = abs(full_fit["mu"] - primary["mu"]) / max(primary["mu"], 1e-9)
            if mu_ratio < 0.30:
                primary = full_fit
        i_mu = primary["i0_global"]
        j    = primary["j_global"]
        sl   = primary["mu"]
        ic   = primary["bgr_ic"]
        r2   = primary["r2"]
        i0   = i_mu  # for lag estimation

        # Refit to get standard error
        _, _, _, _, mu_se = linregress(t[i_mu:j], ln_od[i_mu:j])
        n_exp = j - i_mu

        res["mu"]     = float(sl)
        res["mu_se"]  = float(mu_se)
        res["n_exp"]  = int(n_exp)
        res["td_min"] = float(np.log(2) / sl * 60)
        res["r2"]     = float(r2)
        res["exp_t"]  = t[i_mu:j]
        res["bgr_sl"] = float(sl)
        res["bgr_ic"] = float(ic)

        # ── Quality score (0–100) ────────────────────────────────────────────
        q_r2    = r2 * 40.0
        q_pts   = min(n_exp / 15.0, 1.0) * 30.0
        _resids_final = ln_od[i_mu:j] - (sl * t[i_mu:j] + ic)
        _rstd_final   = float(np.std(_resids_final))
        q_noise = max(0.0, 1.0 - _rstd_final / 0.30) * 20.0
        q_regime = 10.0 if 0.05 <= sl <= 2.0 else 5.0
        res["quality_score"] = int(round(q_r2 + q_pts + q_noise + q_regime))

        # ── Lag: piecewise regression (flat + exponential) breakpoint ────────
        def _detect_lag_piecewise(t_arr, ln_arr, mu_val, min_exp=5, min_lag=2):
            nn = len(t_arr)
            if nn < min_lag + min_exp:
                return None
            best_rss = np.inf
            best_tb  = None
            for bp in range(min_lag, nn - min_exp + 1):
                seg1 = ln_arr[:bp]
                mean1 = np.mean(seg1)
                rss1  = float(np.sum((seg1 - mean1) ** 2))
                seg2_t  = t_arr[bp:]
                seg2_ln = ln_arr[bp:]
                sl2, ic2, _, _, _ = linregress(seg2_t, seg2_ln)
                if sl2 < mu_val * 0.25:
                    continue
                rss2 = float(np.sum((seg2_ln - (sl2 * seg2_t + ic2)) ** 2))
                total_rss = rss1 + rss2
                if total_rss < best_rss:
                    best_rss = total_rss
                    best_tb  = float(t_arr[bp])
            return best_tb

        # For lag, use data up to the end of the primary phase's segment
        seg_end_idx = j
        lag_bp = _detect_lag_piecewise(t[:seg_end_idx], ln_od[:seg_end_idx], sl)

        if lag_bp is not None and lag_bp > t[0]:
            bp_idx = int(np.searchsorted(t[:seg_end_idx], lag_bp))

            flat_od = od[:bp_idx]
            flat_cv = float(np.std(flat_od) / max(np.mean(flat_od), 1e-9))
            is_flat = flat_cv < 0.30 and bp_idx >= 2

            if is_flat:
                res["lag_h"] = float(lag_bp)
            else:
                od_inoc = float(od[0])
                snr_threshold = max(od_inoc * 3.0, 0.3)
                snr_mask = (od[:seg_end_idx] >= snr_threshold)
                snr_idx  = np.where(snr_mask)[0]
                if len(snr_idx) >= MIN_PTS:
                    sl_clean, ic_clean, _, _, _ = linregress(
                        t[snr_idx], ln_od[snr_idx])
                    if sl_clean > 0:
                        od_inoc = max(od_inoc, 1e-9)
                        lag_t = (np.log(od_inoc) - ic_clean) / sl_clean
                        res["lag_h"] = float(max(lag_t, 0.0))
                    else:
                        res["lag_h"] = float(lag_bp)
                else:
                    res["lag_h"] = float(lag_bp)
        elif i0 > 0:
            od_ref = float(np.median(od[:i0]))
            od_ref = max(od_ref, 1e-9)
            lag_t  = (np.log(od_ref) - ic) / sl
            if lag_t > float(t[i0]) and od[i0] >= od[0] * 0.40:
                lag_t = float(t[i0])
            res["lag_h"] = float(max(lag_t, 0.0))
        else:
            od_ref = float(od[0])
            od_ref = max(od_ref, 1e-9)
            lag_t  = (np.log(od_ref) - ic) / sl
            if lag_t < 0:
                lag_t = 0.0
            res["lag_h"] = float(max(lag_t, 0.0))
    else:
        res["mu"] = res["td_min"] = res["r2"] = None
        res["exp_t"] = np.array([])
        res["bgr_sl"] = res["bgr_ic"] = None
        res["lag_h"]  = float(t[0])
        res["quality_score"] = 0
        res["mu_se"] = None
        res["n_exp"] = 0

    return res
