"""Explicit, external-only 2D Meep job specification; does not execute Meep.

2D cylinders are infinite along the out-of-plane direction, and cannot
represent finite 3D spheres or their absolute 3D scattering cross sections.
No automatic network or subprocess invocation is permitted.
"""
from dataclasses import dataclass,asdict
import json
import numpy as np

LIMITATION=("Meep 2D models translationally invariant out-of-plane geometries "
            "(e.g., infinite cylinders), NOT finite 3D spheres. Do not "
            "interpret 2D scattering as 3D sphere Mie cross sections.")

@dataclass
class Meep2DJob:
    wavelength_nm: float
    cylinder_radius_nm: float
    cylinder_index_n: float
    cylinder_index_k: float = 0.
    host_index: float = 1.
    resolution_pixels_per_um: int = 50
    polarization: str = "Ez"
    dimension: int = 2
    user_selected: bool = False

    def validate(self):
        values=(self.wavelength_nm,self.cylinder_radius_nm,self.cylinder_index_n,
                self.cylinder_index_k,self.host_index)
        if not np.all(np.isfinite(values)) or min(values[:3])<=0 or self.host_index<=0 or self.cylinder_index_k<0:
            raise ValueError("Invalid geometry or passive optical constants")
        if self.dimension != 2 or self.polarization not in ("Ez","Hz"):
            raise ValueError("Only explicitly selected 2D Ez/Hz jobs are supported")
        if not isinstance(self.resolution_pixels_per_um,int) or not 10<=self.resolution_pixels_per_um<=500:
            raise ValueError("Resolution must be 10–500 pixels/um")
        if not self.user_selected:
            raise PermissionError("External Meep execution requires explicit user selection")
        return self

    def export_json(self):
        self.validate()
        return json.dumps({"schema":"ocl-meep2d-job-v1","geometry":"infinite_cylinder",
                           "limitations":LIMITATION,"job":asdict(self)},indent=2)
