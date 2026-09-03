"""
Layer and Structure representation for multilayer acoustic simulations.
"""

from dataclasses import dataclass, field
from typing import List, Union, Optional
from pathlib import Path
from .material import Material


@dataclass
class FilmLayer:
    """
    A single film layer in a multilayer structure.

    Attributes:
        material: Material instance or material properties.
        thickness_mean: Mean layer thickness.
        roughness: RMS fluctuation/roughness in layer thickness.
    """
    material: Material
    thickness_mean: float
    roughness: float = 0.0

    def sample_thickness(self, random_normal: float = 0.0) -> float:
        """Sample thickness given a standard normal random value."""
        t = self.thickness_mean + random_normal * self.roughness
        if t <= 0:
            raise ValueError(
                f"Sampled thickness for layer ({self.material.name}) is non-positive: {t}"
            )
        return t


@dataclass
class Structure:
    """
    Multilayer structure consisting of ordered film layers.
    Layer 0 is the surface layer exposed to pump/probe light.
    The last layer is typically the substrate or bottom film.
    """
    layers: List[FilmLayer] = field(default_factory=list)

    def add_layer(self, layer: FilmLayer) -> None:
        """Add a layer to the structure."""
        self.layers.append(layer)

    @property
    def num_layers(self) -> int:
        """Number of layers in the structure."""
        return len(self.layers)

    def sample_structure(self, random_normals: Optional[List[float]] = None) -> List[float]:
        """
        Sample layer thicknesses using random standard normal values.
        If random_normals is None, mean thicknesses are returned.
        """
        if random_normals is None:
            return [layer.thickness_mean for layer in self.layers]
        
        if len(random_normals) != len(self.layers):
            raise ValueError(
                f"Expected {len(self.layers)} random values, got {len(random_normals)}"
            )
            
        return [
            layer.sample_thickness(rnd)
            for layer, rnd in zip(self.layers, random_normals)
        ]
