#!/usr/bin/env python3
"""Freeze prospective flight outputs and reference-independent transport checks."""
import argparse
import json
from pathlib import Path
import numpy as np
from prepare_prospective_flight import ROOT, OUT
from research_flight_conditional import load, write, digest
from research_conditional_mc import stamp
import flight_conditional_likelihood as conditional
import joint_adc_likelihood as joint
import constrained_joint_adc as constrained
import free_timing_transport as center
import continuous_rich_transport as forward
import rich_transport as material

JOINT = joint.OUT / 'flight'
FAMILIES = [('native_gaussianR', 't_native'), ('newcalib_gaussianR', 't'), ('native_studentR5', 't_native')]
BRANCHES = ['ols', 'time', 'joint', 'speed_null', 'speed_gain_null']
MODES = ['fit', 'unit']
SCALES = [('s1', 1.), ('s0p5', .5), ('s2', 2.)]
SOURCES = ['predict_prospective_flight.py', 'joint_adc_likelihood.py', 'constrained_joint_adc.py',
           'flight_conditional_likelihood.py', 'free_timing_transport.py', 'flight_mean_loss.py',
           'research_conditional_output.py', 'continuous_rich_transport.py', 'rich_transport.py']


def plan():
    path = OUT / 'prediction_protocol.json'
    if path.exists():
        raise FileExistsError('Prospective output protocol already exists')
    model_hash = {name: digest(JOINT / name / 'freeze.json') for name, _ in FAMILIES}
    primary_baseline = 'newcalib_gaussianR_ols_fit_RICH_s1_u'
    write(path, dict(created_utc=stamp(), families=FAMILIES, branches=BRANCHES, material_modes=MODES, rich_material_scales=SCALES,
        source_sha256={name: digest(ROOT / 'analysis' / name) for name in SOURCES}, frozen_model_sha256=model_hash,
        normalization_sha256=digest(JOINT / 'protocol.json'), material_map_sha256=digest(material.MAP),
        measured_selection_sha256=digest(OUT / 'native_domain_protocol.json'),
        candidates_sha256=digest(OUT / 'candidates.npz'),
        primary_comparisons=[['newcalib_gaussianR_joint_fit_RICH_s1_u', primary_baseline],
                             ['newcalib_gaussianR_speed_gain_null_fit_RICH_s1_u', primary_baseline]],
        primary_endpoint='W68 of reconstructed beta at the NaF emission plane minus the independently reconstructed RICH beta, on the identical selected reference events',
        primary_scope='NaF strict existing reference quality, measured native total Tracker chi-square <=10 in addition to the previously declared X/Y cuts; no reconstructed beta, mass, or TOF/RICH agreement cut',
        inference='Original separately flight-trained 128-node joint models; no retraining or response recentering on the new measured events',
        controls='Every family includes OLS, marginal timing GLS, joint time/ADC, ADC-speed-null and ADC-speed-plus-common-gain-null output. Native BetaH and constant beta=1 are retained as untransported controls, not called L5 estimates.',
        mean_loss='Each candidate and its matching OLS/time baseline use the same fitted two-group first-order L5 inverse. Unit-material L5 inverse and half/double forward material are declared sensitivities, not a posterior branch choice.',
        rich_transport='Existing elemental g4mscan column map and Tracker Hermite geometry; midpoint continuous-mass Bethe mean propagation from L5 to the geometrically predicted NaF emission z=-75.552 cm',
        beta_above_one='u<=1 is retained unchanged by the forward material map, with zero mean-loss extension, not clipped or removed. No isotope/source posterior enters inference.',
        propagation_failure='Preserve invalid-output flags and all input identities; do not exclude an event based on which estimator propagated. A branch with any failed output in its reference cohort fails complete-coverage acceptance.',
        numerical_tests='Five one-way transport tests plus five constrained-noise tests; fixed first 128 measured events for 1/2/4-cell propagation convergence check before opening RICH values',
        response='Separate exact time-only response of constrained/OLS/time outputs from the joint model combined time+ADC response; verify the existing material-aware speed injection for center inference',
        evaluation_bootstrap=dict(replicates=4096, unit='entire new evaluation run', seed=1003202610,
                                  quantiles=[2.5, 50, 97.5], two_primary_adjusted_quantiles=[1.25, 50, 98.75]),
        recalibration_bootstrap='Use all 128 previously saved flight NewCalib/Gaussian run refits and their own feature normalization; 16 paired new-evaluation-run draws per calibration. Original 128-node and 256-node refit point difference must remain explicit.',
        bins=dict(measured_R=[2.5, 3, 3.5, 4.5, 6], reference_beta=[.78, .82, .86, .9, .94, .98, 1.02],
                  year=[2011, 2014, 2017, 2020, 2023, 2027]),
        new_RICH_values_read=False, MC_coefficients_used=False))
    print('PROSPECTIVE_OUTPUT_PROTOCOL', flush=True)


def check():
    p = json.loads((OUT / 'prediction_protocol.json').read_text())
    for name, fingerprint in p['source_sha256'].items():
        if digest(ROOT / 'analysis' / name) != fingerprint:
            raise ValueError('Prediction source changed ' + name)
    for name, fingerprint in p['frozen_model_sha256'].items():
        if digest(JOINT / name / 'freeze.json') != fingerprint:
            raise ValueError('Flight model changed')
    if digest(JOINT / 'protocol.json') != p['normalization_sha256']:
        raise ValueError('Flight normalization changed')
    if digest(OUT / 'native_domain_protocol.json') != p['measured_selection_sha256']:
        raise ValueError('Measured input domain changed')
    return p


def geometry():
    p = check()
    if (OUT / 'rich_geometry.npz').exists():
        raise FileExistsError('Prospective reference-plane geometry already exists')
    if digest(material.MAP) != p['material_map_sha256']:
        raise ValueError('Elemental material map changed')
    a = load(OUT / 'measured.npz')
    mapper = material.MaterialMap()
    result = mapper.build_columns(a['tk_pos'], a['tk_dir'])
    if not result['valid'].all():
        raise ValueError('Invalid measured RICH-path geometry: %d' % (~result['valid']).sum())
    if not result['NaF_predicted'].all():
        raise ValueError('Central-trajectory selection and NaF geometry are inconsistent')
    np.savez_compressed(OUT / 'rich_geometry.npz', run=a['run'], event=a['event'], **result)
    write(OUT / 'rich_geometry_verification.json', dict(complete=True, created_utc=stamp(), events=len(a['run']), valid=int(result['valid'].sum()),
        measured_sha256=digest(OUT / 'measured.npz'), geometry_sha256=digest(OUT / 'rich_geometry.npz'),
        target_z_values=np.unique(result['target_z']), total_column_percentiles=np.percentile(result['column'].sum((1, 2)), [0, 5, 50, 95, 100]),
        material_map_sha256=p['material_map_sha256'], new_RICH_values_read=False))
    print('PROSPECTIVE_RICH_GEOMETRY', len(a['run']), flush=True)


def infer_family(theta, data, geometry, feature_normalization, field):
    n = len(data['run']); x = geometry['x']
    features, _ = joint.response_features(data, np.ones(n, bool), feature_normalization)
    y = joint.observation(data, data['adc'], field)
    z = joint.reconstruct(theta, x, y, features)
    xc = x - x.mean(1)[:, None]; weights = xc / np.sum(xc * xc, axis=1)[:, None]
    ols = np.sum(weights * data[field], axis=1)
    outputs = dict(ols=(ols, weights), time=(z['time_u'], z['time_weights']),
                   joint=(z['joint_u'], z['weights'][:, :3].dot(conditional.D)))
    responses = dict(joint_time_gain_percentiles=np.percentile(z['time_injection_gain'], [0, 5, 50, 95, 100]),
                     joint_ADC_speed_gain_percentiles=np.percentile(z['adc_injection_gain'], [0, 5, 50, 95, 100]))
    for variant in ['speed_null', 'speed_gain_null']:
        r = constrained.reconstruct(theta, x, y, features, variant)
        outputs[variant] = (r['u'], r['time_weights'])
        responses[variant] = dict(time_gain_max_error=float(np.max(abs(np.sum(r['time_weights'] * x, axis=1) - 1))),
            common_time_max_error=float(np.max(abs(r['time_weights'].sum(1)))),
            constraint_max_error=float(np.max(abs(r['response'] - np.r_[1., np.zeros(r['response'].shape[1] - 1)]))))
        if max(responses[variant].values()) > 2e-10:
            raise ValueError('Prospective constrained timing response failed')
    return outputs, responses


def predict():
    p = check()
    if (OUT / 'predictions_before_reference.npz').exists() or (OUT / 'output_freeze.json').exists():
        raise FileExistsError('Prospective output snapshot already exists')
    data = load(OUT / 'measured.npz'); g = load(OUT / 'moments.npz'); rg = load(OUT / 'rich_geometry.npz')
    for key in ['run', 'event']:
        if not np.array_equal(data[key], g[key]) or not np.array_equal(data[key], rg[key]):
            raise ValueError('Prospective geometry identity mismatch')
    if not g['valid'].all() or not rg['valid'].all():
        raise ValueError('Resolve invalid measured geometry before prediction freeze')
    cp = json.loads((JOINT / 'protocol.json').read_text()); norm = cp['feature_normalization']
    prediction = {key: data[key] for key in ['run', 'event', 'source', 'entry', 'R']}
    prediction['native_BetaH_u'] = 1 / data['tof_betah']; prediction['constant_one_u'] = np.ones(len(data['run']))
    diagnostics = {}; convergence = {}; invalid = {}; ncheck = min(128, len(data['run'])); ix = np.arange(ncheck)
    small_geom = forward.subset_geometry(rg, ix)
    for family, field in FAMILIES:
        theta = np.array(json.loads((JOINT / family / 'freeze.json').read_text())['model']['theta'])
        branches, response = infer_family(theta, data, g, norm, field); diagnostics[family] = response
        for branch, (slope, weight) in branches.items():
            if not np.isfinite(slope).all() or np.any(slope == 0):
                raise ValueError('Nonfinite or zero measured slope; do not remove output-selected events')
            prediction[family + '_' + branch + '_slope_u'] = slope
            for mode in MODES:
                kappa = theta[:2] if mode == 'fit' else np.ones(2)
                center_result = center.solve_center(slope, data['R'], weight, g['moment'], kappa)
                if not center_result['converged'].all():
                    raise ValueError('Prospective mean-loss inverse failed; preserve the measured population')
                u = center_result['u']; label = family + '_' + branch + '_' + mode
                prediction[label + '_L5_u'] = u
                for scale_name, scale in SCALES:
                    result = forward.propagate(u, data['R'], rg, material_scale=scale)
                    key = label + '_RICH_' + scale_name + '_u'
                    prediction[key] = result['u']; prediction[key + '_valid'] = result['valid']
                    invalid[key] = int((~result['valid']).sum())
                    if mode == 'fit' and scale_name == 's1' and family == 'newcalib_gaussianR':
                        one = result['u'][ix]
                        two = forward.propagate(u[ix], data['R'][ix], forward.subdivide(small_geom, 2))
                        four = forward.propagate(u[ix], data['R'][ix], forward.subdivide(small_geom, 4))
                        good = result['valid'][ix] & two['valid'] & four['valid']
                        convergence[branch] = dict(events=ncheck, valid=int(good.sum()),
                            max_u_difference_1_to_2=float(np.max(abs(one[good] - two['u'][good]))),
                            max_u_difference_2_to_4=float(np.max(abs(two['u'][good] - four['u'][good]))),
                            median_material_shift=float(np.median(result['u'][result['valid']] - u[result['valid']])))
                        if convergence[branch]['max_u_difference_1_to_2'] > 1e-5:
                            raise ValueError('Reference-plane material subdivision check exceeds tolerance')
        print('PROSPECTIVE_FAMILY_PREDICTED', family, len(data['run']), flush=True)
    np.savez_compressed(OUT / 'predictions_before_reference.npz', **prediction)
    write(OUT / 'output_freeze.json', dict(complete=True, created_utc=stamp(), events=len(data['run']), runs=len(np.unique(data['run'])),
        prediction_sha256=digest(OUT / 'predictions_before_reference.npz'), protocol_sha256=digest(OUT / 'prediction_protocol.json'),
        measured_sha256=digest(OUT / 'measured.npz'), moments_sha256=digest(OUT / 'moments.npz'),
        rich_geometry_sha256=digest(OUT / 'rich_geometry.npz'), frozen_model_sha256=p['frozen_model_sha256'],
        response_diagnostics=diagnostics, numerical_propagation_checks=convergence, invalid_outputs=invalid,
        new_RICH_values_read=False, no_output_based_event_removal=True, MC_coefficients_used=False))
    print('PROSPECTIVE_OUTPUTS_FROZEN', len(data['run']), 'propagation_failures', sum(invalid.values()), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['plan', 'geometry', 'predict']); args = parser.parse_args()
    globals()[args.stage]()
