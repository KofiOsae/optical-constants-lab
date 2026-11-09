"""Validated 2D/3D Meep geometry job and self-contained Colab notebook export.

First implementation deliberately supports nondispersive, lossless dielectric
objects only. Broadband lossy dispersive materials require fitted causal
susceptibilities and must not be silently approximated by complex epsilon.
"""
from dataclasses import dataclass,asdict
import json,math

@dataclass
class MeepJob:
    dimension:int=2
    geometry:str="cylinder"
    wavelength_nm:float=550.
    radius_nm:float=50.
    particle_n:float=2.
    host_n:float=1.
    resolution_pixels_per_um:int=50
    padding_um:float=1.
    pml_um:float=.5
    polarization:str="Ez"

    def validate(self):
        if self.dimension not in (2,3):raise ValueError("Dimension must be 2 or 3")
        if self.geometry != ("cylinder" if self.dimension==2 else "sphere"):
            raise ValueError("2D supports infinite cylinders; 3D supports finite spheres")
        vals=(self.wavelength_nm,self.radius_nm,self.particle_n,self.host_n,self.padding_um,self.pml_um)
        if any(not math.isfinite(v) or v<=0 for v in vals):raise ValueError("Positive finite geometry and material parameters required")
        if self.polarization not in ("Ez","Hz") and self.dimension==2:raise ValueError("2D polarization must be Ez or Hz")
        if not isinstance(self.resolution_pixels_per_um,int) or not 10<=self.resolution_pixels_per_um<=200:
            raise ValueError("Resolution must be 10–200 pixels/um")
        if self.radius_nm/1000 < 4/self.resolution_pixels_per_um:
            raise ValueError("Radius resolved by fewer than 4 pixels; increase resolution")
        side=2*(self.radius_nm/1000+self.padding_um+self.pml_um)
        cells=(side*self.resolution_pixels_per_um)**self.dimension
        if cells>2e7:raise ValueError("Estimated grid exceeds research preview safety cap (20 million cells)")
        return self

    def job_dict(self):
        self.validate()
        return {"schema":"ocl-meep-job-v1","dimension":self.dimension,
                "geometry":self.geometry,"units":"um","job":asdict(self),
                "limitations":"lossless nondispersive dielectric, external Colab execution; no optical cross-section extraction yet"}

def make_notebook(job):
    job.validate()
    config=json.dumps(job.job_dict(),indent=2)
    cells=[]
    def md(src):cells.append({"cell_type":"markdown","metadata":{},"source":src.splitlines(True)})
    def code(src):cells.append({"cell_type":"code","execution_count":None,"metadata":{},"outputs":[],"source":src.splitlines(True)})
    md("# Optical Constants Laboratory · Meep research simulation\n"
       "This notebook is preconfigured. **Review the inputs and authorize Colab before running.** "
       "2D is an infinite cylinder; 3D is a finite sphere. Results are not equivalent. "
       "This prototype produces fields, not validated scattering cross sections.")
    md("## Runtime setup · Colab-managed conda-forge PyMeep\\n"
       "The next cell attempts installation using conda-forge, not the unrelated PyPI package named meep. "
       "Colab may restart its kernel during conda setup. If it does, select Runtime → Run all again. "
       "The runtime must be authorized by the user. Installation is network- and Colab-version-dependent; "
       "a fresh-runtime success has not been independently verified.")
    code("""import importlib.util, subprocess, sys
if importlib.util.find_spec('meep') is None:
    try:
        import condacolab
    except ImportError:
        subprocess.check_call([sys.executable,'-m','pip','install','-q','condacolab'])
        import condacolab
    condacolab.install()
    subprocess.check_call(['mamba','install','-y','-c','conda-forge','pymeep','numpy','matplotlib'])
import meep as mp
import numpy as np
import matplotlib.pyplot as plt
print('PyMeep import successful:',getattr(mp,'__version__','unknown'))
""")
    code("import json\nconfig = json.loads("+repr(config)+")\nprint('Loaded job:',config['schema'],config['dimension'],'D')")
    code("""j=config['job']
r=j['radius_nm']/1000.
side=2*(r+j['padding_um']+j['pml_um'])
dim=j['dimension']
cell=mp.Vector3(side,side,0 if dim==2 else side)
material=mp.Medium(index=j['particle_n'])
background=mp.Medium(index=j['host_n'])
geometry=[mp.Cylinder(radius=r,height=mp.inf,material=material)] if dim==2 else [mp.Sphere(radius=r,material=material)]
component=(mp.Ez if j['polarization']=='Ez' else mp.Hz) if dim==2 else mp.Ez
source=mp.Source(mp.GaussianSource(frequency=1000./j['wavelength_nm'],fwidth=.15*1000./j['wavelength_nm']),
    component=component,center=mp.Vector3(-side/2+j['pml_um']+.25,0,0),
    size=mp.Vector3(0,side-2*j['pml_um'],0 if dim==2 else side-2*j['pml_um']))
sim=mp.Simulation(cell_size=cell,boundary_layers=[mp.PML(j['pml_um'])],
    geometry=geometry,default_material=background,sources=[source],
    resolution=j['resolution_pixels_per_um'],dimensions=dim)
sim.run(until_after_sources=100)
field=np.asarray(sim.get_array(center=mp.Vector3(),size=mp.Vector3(side-2*j['pml_um'],side-2*j['pml_um'],0),component=component))
plt.figure(figsize=(6,5));plt.imshow(field.T,origin='lower',cmap='RdBu');plt.colorbar(label=str(component));plt.title(f'{dim}D field slice');plt.show()
""")
    md("## Assisted numerical assessment\\nRepeat the simulation with higher resolution, larger padding and PML. "
       "A field snapshot is not a scattering cross section or a converged quantitative observable.")
    code("""print('Meep version:',getattr(mp,'__version__','unknown'))
print('Geometry:',j['geometry'],'dimension:',dim,'resolution:',j['resolution_pixels_per_um'])
print('Estimated grid cells:',round((side*j['resolution_pixels_per_um'])**dim))
print('Field shape:',field.shape,'finite:',bool(np.all(np.isfinite(field))))
if not np.all(np.isfinite(field)):
    raise ValueError('Nonfinite simulated fields; check PML and numerical stability')
print('Review: resolution convergence, PML reflection, background reference and incident normalization.')
""")
    code("""import csv
from google.colab import files
name='ocl_meep_field_slice.csv'
np.savetxt(name,field,delimiter=',')
metadata={'schema':'ocl-meep-field-v1','dimension':dim,'geometry':j['geometry'],
          'wavelength_nm':j['wavelength_nm'],'field_component':str(component),
          'resolution_pixels_per_um':j['resolution_pixels_per_um'],
          'data_file':name,'note':'Field snapshot only; not R/T or normalized scattering'}
with open('ocl_meep_field_metadata.json','w') as f:json.dump(metadata,f,indent=2)
files.download(name);files.download('ocl_meep_field_metadata.json')
""")
    return {"nbformat":4,"nbformat_minor":5,"metadata":{"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"},
            "colab":{"name":"ocl_meep_research.ipynb"}},"cells":cells}
