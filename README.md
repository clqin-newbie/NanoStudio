# [NanoStudio: A Python Package for Automated Construction of Twisted Two-Dimensional Materials and Multi-Walled Nanotubes]()

Here, we present NanoStudio, an open-source Python package for the automated construction of low-dimensional nanostructures based on honeycomb lattices. The package provides a unified implementation for generating graphene with arbitrary chirality, multilayer graphene with user-defined stacking sequences, twisted graphene with arbitrary commensurate twist angles, quasi-one-dimensional nanoribbons, and multi-walled nanotubes. 

![alt text](image.png)

## Authors

This package was primarily written by Chenglong Qin (clqin@xhu.edu.cn).

## License

The model is released under the MIT License.

Parameter Documentation for Nanostructure Generation CodeThis object-oriented code generates graphene, twisted bilayer graphene (TBG), and carbon nanotubes (CNTs) based on the ASE library. Below is a full English summary of classes, parameters, methods and usage examples.1. Base Class: NanoAbstract base class providing common parameters, coordinate transformation, mathematical tools and basic structure functions for all carbon nanostructures.1.1 Constructor __init__()




















































表格ParameterTypeDescriptionUnit & NotenintChiral index nStandard index for carbon nanostructuresmintChiral index mStandard index for carbon nanostructuresperiodlistSupercell repetition [x, y, z]Lattice replication along three axespbclist[bool]Periodic boundary conditions [x, y, z]True: enable PBC; False: add vacuum layerbond_lengthfloatC-C bond lengthUnit: Å, default = 1.42 Åatom_typelistElement symbol of atomse.g. ['C'], ['C', 'B']vacuumfloatVacuum layer thicknessUnit: Å, applied to non-periodic directions1.2 General Methods
rotation_a1a2(theta): 2D rotation matrix under oblique lattice basis (a1, a2). Input theta in radians.
a1a2_to_ij(theta): Transformation matrix from oblique coordinates to Cartesian coordinates.
rotation_ij(theta): Standard 2D rotation matrix in Cartesian system. Input theta in radians.
extended_gcd(a, b): Extended Euclidean algorithm. Solve \(ax+by=\gcd(a,b)\), returns (gcd, x, y).
_get_structures(): Generate ASE Atoms object with lattice, atomic positions, PBC and vacuum settings.
_get_properties(): Calculate structural properties including chiral angle, bond length and PBC.
show_properties(): Print structural parameters in console.
2. Subclass: Graphene(Nano)For constructing single-layer and multi-layer stacked graphene (AA / AB / arbitrary stacking).2.1 Constructor __init__()
































































表格ParameterTypeDescriptionDefault ValuenintChiral index n—mintChiral index m—pintComponent p of translation vector T—qintComponent q of translation vector T—xy_periodlistSupercell replication along x & y [nx, ny][1, 1]xy_pbclist[bool]PBC for x & y directions[True, True]bond_lengthfloatC-C bond length1.42 Åatom_typelistAtomic element['C']vacuumfloatVacuum layer thickness7.5 Å
Internal logic: Z direction is automatically set to non-periodic with no replication.
2.2 Main Method get_structures()




























表格ParameterTypeDescriptionRestrictionlayer_numintNumber of graphene layersMust match length of translationstranslationslist[tuple]In-plane offset for each layer (fractional coordinates)Format: [(x1,y1), (x2,y2), ...]layer_spacingfloatVertical interlayer distanceUnit: Å, typical value = 3.4 Å2.3 Private Methods
_build_coords(): Generate fractional coordinates for single-layer graphene (A / B sublattice).
_get_properties(): Calculate lattice matrix, total atom number, chiral angle and translation vectors.
3. Subclass: TwistGraphene(Nano)For constructing twisted multi-layer graphene (TBG) with rotational misalignment between layers.3.1 Constructor __init__()




















































表格ParameterTypeDescriptionDefault ValuenintChiral index n—mintChiral index m—xy_periodtupleSupercell replication (nx, ny)(1, 1)xy_pbclist[bool]PBC for x & y directions[True, True]bond_lengthfloatC-C bond length1.42 Åatom_typelistAtomic element['C']vacuumfloatVacuum layer thickness7.5 Å3.2 Main Method get_structures()


































表格ParameterTypeDescriptionRestrictionlayer_numintTotal number of layersLength must match translations & layer_sequencetranslationslist[tuple]In-plane offset for each layerFractional coordinateslayer_sequencelist[int]Layer type flag0: original layer; 1: rotated layerlayer_spacingfloatVertical interlayer distanceUnit: Å3.3 Notes
Twist angle is automatically calculated from (n, m).
Layer rotation is realized via Cartesian rotation matrix.
_get_single_layer(): Generate coordinates for original or rotated graphene monolayer.
4. Subclass: Nanotube(Nano)For constructing single-walled and multi-walled carbon nanotubes (CNTs); supports multiple tubes in one model.4.1 Constructor __init__()




















































表格ParameterTypeDescriptionDefault Valuenlist[int]List of chiral index n for each tubee.g. [15, 6] for two tubesmlist[int]List of chiral index m for each tubeSame length as nz_periodlist[int]Axial replication along tube axis (z)One value per tubez_pbcboolPBC along z-axis (tube axis)Truebond_lengthfloatC-C bond length1.42 Åatom_typelistAtomic element['C']vacuumfloatRadial vacuum thickness7.5 Å
Rule: PBC is disabled for radial x & y directions by default.
4.2 Main Method get_structures()
















表格ParameterTypeDescriptiontranslationlist[tuple]In-plane offset for each nanotubeLength equals total tube number4.3 Calculated Structural Properties
chirality (n, m): Chiral indices
diameter: Tube diameter (Å)
length: Axial period length (Å)
chiral_theta: Chiral angle (°)
total_atom_num: Total number of atoms
4.4 Coordinate Logic
2D graphene plane is rolled into cylindrical coordinates to form CNT geometry.
Atomic positions for A / B sublattices are calculated via rotation angle and axial shift.
5. Default Global Parameters





























表格Physical QuantityDefault ValueRemarksC-C bond length1.42 ÅStandard value for carbon materialsVacuum layer7.5 ÅUsed for all non-periodic directionsGraphene interlayer distance3.4 ÅTypical van der Waals gap for graphiteDefault atom type['C']Pure carbon structure6. Code Usage Examples6.1 Multilayer Graphenepython运行model = Graphene(4, 2, -2, 6, xy_period=[1, 1], xy_pbc=[True, True])
struct = model.get_structures(layer_num=3, translations=[(0, 0), (1/3, 1/3), (-1/3, 2/3)])

Chirality: \(n=4, m=2\); Translation vector: \(p=-2, q=6\)
3 layers with different in-plane shifts for varied stacking configurations
6.2 Twisted Graphenepython运行model = TwistGraphene(4, 2, xy_period=[2, 2], xy_pbc=[True, True])
struct = model.get_structures(layer_num=3, translations=[(0, 0), (0, 0), (0, 0)], layer_sequence=[1, 0, 1])

Chirality: \(n=4, m=2\); 2×2 supercell in x-y plane
3 layers: rotated layer → original layer → rotated layer
6.3 Multi-Wall Carbon Nanotubespython运行model = Nanotube([15, 6], [15, 6], z_period=[5, 5], atom_type=['C'], z_pbc=False)
struct = model.get_structures(translation=[(0, 0), (0, 0)])

Two armchair CNTs: \((15,15)\) and \((6,6)\)
5 times axial replication; PBC disabled along tube axis; no in-plane offset