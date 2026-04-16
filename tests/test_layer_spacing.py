import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from twist_generator.build import build_bilayer_aligned
from twist_generator.structure import atomic_structure, atom_unit
from twist_generator.tools import combine_lattice


def make_lattice(c, z_cartesian, element='X', xy=None):
    lattice = atomic_structure(mode='empty')
    lattice.title = 'test'
    lattice.frac = 1.0
    lattice.basis = [[1., 0., 0.], [0., 1., 0.], [0., 0., c]]
    lattice.elements = [element]
    lattice.natoms = [len(z_cartesian)]
    lattice.nions = len(z_cartesian)
    lattice.selective_dynamics = 0
    lattice.tag = 'Direct'
    if xy is None:
        xy = [(0.1 + 0.1*i, 0.2 + 0.1*i) for i in range(len(z_cartesian))]
    lattice.atom = [
        atom_unit([x, y, z/c], element)
        for (x, y), z in zip(xy, z_cartesian)
    ]
    return lattice


def z_positions_cartesian(lattice):
    c = lattice.cell_length()[-1]
    return [atom.position[2]*c for atom in lattice.atom]


class LayerSpacingTests(unittest.TestCase):
    def assert_close(self, actual, expected, places=8):
        self.assertAlmostEqual(actual, expected, places=places)

    def test_toy_slab_keeps_surface_gap_and_original_external_vacuum(self):
        lattice = make_lattice(10., [4.5, 5.5])

        bilayer = combine_lattice(lattice, lattice, 4.)
        z = z_positions_cartesian(bilayer)
        layer_a = z[:2]
        layer_b = z[2:]

        self.assert_close(bilayer.cell_length()[-1], 15.)
        self.assert_close(min(layer_b)-max(layer_a), 4.)
        self.assert_close(min(z)+(bilayer.cell_length()[-1]-max(z)), 9.)

    def test_heterobilayer_uses_larger_input_vacuum(self):
        lattice_a = make_lattice(20., [9., 11.], element='A')
        lattice_b = make_lattice(12., [4.5, 7.5], element='B')

        bilayer = combine_lattice(lattice_a, lattice_b, 3.3)
        z = z_positions_cartesian(bilayer)
        layer_a = z[:2]
        layer_b = z[2:]

        self.assert_close(bilayer.cell_length()[-1], 26.3)
        self.assert_close(min(layer_b)-max(layer_a), 3.3)
        self.assert_close(min(z)+(bilayer.cell_length()[-1]-max(z)), 18.)

    def test_wrapped_z_slab_uses_largest_periodic_gap(self):
        lattice = make_lattice(10., [9.5, 0.5])

        bilayer = combine_lattice(lattice, lattice, 4.)
        z = z_positions_cartesian(bilayer)
        layer_a = z[:2]
        layer_b = z[2:]

        self.assert_close(bilayer.cell_length()[-1], 15.)
        self.assert_close(max(layer_a)-min(layer_a), 1.)
        self.assert_close(min(layer_b)-max(layer_a), 4.)

    def test_vacuum_override_sets_final_external_vacuum(self):
        lattice = make_lattice(10., [4.5, 5.5])

        bilayer = combine_lattice(lattice, lattice, 4., vacuum=15.)
        z = z_positions_cartesian(bilayer)

        self.assert_close(bilayer.cell_length()[-1], 21.)
        self.assert_close(min(z)+(bilayer.cell_length()[-1]-max(z)), 15.)

    def test_aligned_build_applies_xy_translation_to_second_layer(self):
        lattice_a = make_lattice(10., [5.], element='A', xy=[(0.1, 0.2)])
        lattice_b = make_lattice(10., [5.], element='B', xy=[(0.1, 0.2)])

        bilayer = build_bilayer_aligned(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            4.,
            lattice_a,
            lattice_b,
            translation=[0.25, 0.1, 0.],
        )

        z = z_positions_cartesian(bilayer)
        self.assert_close(bilayer.atom[1].position[0], 0.35)
        self.assert_close(bilayer.atom[1].position[1], 0.3)
        self.assert_close(z[1]-z[0], 4.)

    def test_invalid_distance_and_vacuum_raise(self):
        lattice = make_lattice(10., [4.5, 5.5])

        with self.assertRaises(ValueError):
            combine_lattice(lattice, lattice, -1.)
        with self.assertRaises(ValueError):
            combine_lattice(lattice, lattice, 4., vacuum=-1.)

    def test_tilted_c_axis_raises(self):
        lattice = make_lattice(10., [4.5, 5.5])
        tilted = make_lattice(10., [4.5, 5.5])
        tilted.basis[2] = [0.1, 0., 10.]

        with self.assertRaises(ValueError):
            combine_lattice(lattice, tilted, 4.)

    def test_bulk_like_distribution_without_unique_vacuum_raises(self):
        lattice = make_lattice(10., [0., 2.5, 5., 7.5])

        with self.assertRaises(ValueError):
            combine_lattice(lattice, lattice, 4.)


if __name__ == '__main__':
    unittest.main()
