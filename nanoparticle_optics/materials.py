"""Metadata-only research material records; no production database changes.

All wavelengths in nm. Provenance and license must be explicit for exports.
"""
from dataclasses import dataclass, field, asdict
import json
import numpy as np


@dataclass
class MaterialRecord:
    material_id: str
    name: str
    wavelength_nm: list
    n: list
    k: list
    composition: dict = field(default_factory=dict)
    dopants: dict = field(default_factory=dict)
    particle_radius_nm: float | None = None
    particle_radius_std_nm: float | None = None
    synthesis: dict = field(default_factory=dict)
    measurement: dict = field(default_factory=dict)
    provenance: str = ""
    data_license: str = ""
    notes: str = ""

    def validate(self):
        if not self.material_id.strip() or not self.name.strip():
            raise ValueError("Material ID and name required")
        wl=np.asarray(self.wavelength_nm, float)
        n=np.asarray(self.n,float)
        k=np.asarray(self.k,float)
        if (wl.ndim != 1 or len(wl)<2 or n.shape != wl.shape or k.shape != wl.shape
            or not np.all(np.isfinite(wl)) or not np.all(np.diff(wl)>0)
            or np.any(wl<=0) or not np.all(np.isfinite(n)) or np.any(n<0)
            or not np.all(np.isfinite(k)) or np.any(k<0)):
            raise ValueError("Invalid passive n,k spectrum or wavelength grid")
        for value in (self.particle_radius_nm, self.particle_radius_std_nm):
            if value is not None and (not np.isfinite(value) or value < 0):
                raise ValueError("Particle radius metadata must be nonnegative")
        if self.particle_radius_nm is not None and self.particle_radius_nm <= 0:
            raise ValueError("Particle radius must be positive")
        if not self.provenance.strip() or not self.data_license.strip():
            raise ValueError("Provenance and data license required")
        return self

    def to_json(self):
        self.validate()
        return json.dumps(asdict(self),indent=2,allow_nan=False)

    @classmethod
    def from_json(cls, text):
        record=cls(**json.loads(text))
        return record.validate()

    def complex_index(self):
        self.validate()
        return np.asarray(self.n)+1j*np.asarray(self.k)

    def permittivity(self):
        return self.complex_index()**2
