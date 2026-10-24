"""Compare the opt-in local quadrature transfer on issue #280's edge over box.

Run with native IPC Python bindings, numpy and matplotlib:
    python esp_local_transfer.py --output results
Forces are analytic derivatives; the reported FD check rebuilds contacts.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import ipctk


def scene():
    rest = np.array([[.035, -.016, .0006], [.065, .017, .0008],
                     [-.05, -.04, -.1], [.05, -.04, -.1], [.05, .04, -.1], [-.05, .04, -.1],
                     [-.05, -.04, 0], [.05, -.04, 0], [.05, .04, 0], [-.05, .04, 0]])
    faces = np.array([[0,3,2], [0,2,1], [4,5,6], [4,6,7], [0,1,5], [0,5,4],
                      [3,7,6], [3,6,2], [0,4,7], [0,7,3], [1,2,6], [1,6,5]]) + 2
    edges = np.unique(np.sort(faces[:, ((0,1), (1,2), (2,0))].reshape(-1,2), axis=1), axis=0)
    mesh = ipctk.CollisionMesh(rest, np.vstack(([[0,1]], edges)), faces)
    mesh.can_collide = ipctk.make_vertex_patches_filter(np.r_[np.zeros(2, dtype=int), np.ones(8, dtype=int)])
    mesh.init_adjacencies()
    return mesh, rest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('local-transfer-results'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    mesh, rest = scene()
    params, potentials = [], []
    for enabled in (False, True):
        p = ipctk.ESPParameters(.03, .05, area_weights=True)
        p.local_quadrature_transfer = enabled
        params.append(p)
        potentials.append(ipctk.ESPPotential(p, use_near_far=True))
    contacts = ipctk.ESPCollisions()
    def evaluate(x, mode):
        X = rest.copy()
        X[0,0] = x
        contacts.build(mesh, X, params[mode])
        e = potentials[mode](contacts, mesh, X)
        f = -potentials[mode].gradient(contacts, mesh, X).reshape(-1,3)
        return np.r_[e, f[0]]
    xs = np.linspace(.04970, .050015, 631)
    data = np.array([np.r_[x*1000, evaluate(x,0), evaluate(x,1)] for x in xs])
    np.savetxt(args.output/'comparison.csv', data, delimiter=',', comments='',
               header='Ax_mm,old_energy,old_Fx,old_Fy,old_Fz,new_energy,new_Fx,new_Fy,new_Fz')
    summary = {'samples': len(xs), 'range_m': [float(xs[0]),float(xs[-1])], 'mollifier': 'unchanged', 'modes': {}}
    fig, axes = plt.subplots(3, 1, figsize=(9,8), sharex=True, layout='constrained')
    for mode, label, color, start in [(0,'Original','#b84358',1),(1,'Local transfer','#247f92',5)]:
        peak = np.argmax(np.abs(data[:,start+1]))
        h = 1e-9
        checked = []
        for x in [.04993,.04995,.0499615,.049984,.049995,.05001]:
            v = evaluate(x,mode)
            fd = -(evaluate(x+h,mode)[0]-evaluate(x-h,mode)[0])/(2*h)
            np.testing.assert_allclose(fd,v[1],rtol=2e-5,atol=2e-5)
            checked.append({'x_m':x,'analytic_Fx':float(v[1]),'FD_Fx':float(fd)})
        summary['modes'][label] = {'peak_abs_Fx':float(abs(data[peak,start+1])),
            'peak_x_mm':float(data[peak,0]),'at_original_peak':evaluate(.0499615,mode).tolist(), 'force_checks':checked}
        for ax,col in zip(axes,[start+1,start+3,start]):
            ax.plot(data[:,0],data[:,col],label=label,color=color)
    summary['peak_reduction_percent'] = 100*(1-summary['modes']['Local transfer']['peak_abs_Fx']/summary['modes']['Original']['peak_abs_Fx'])
    for ax,ylabel in zip(axes,['Endpoint A horizontal force','Endpoint A vertical force','ESP energy']):
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.2)
        ax.legend(frameon=False)
        ax.ticklabel_format(axis='x',style='plain',useOffset=False)
    axes[-1].set_xlabel('Endpoint A.x (mm)')
    fig.suptitle('Issue #280 · edge above box · current mollifier in both modes')
    fig.savefig(args.output/'comparison.png',dpi=180)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__ == '__main__':
    main()
