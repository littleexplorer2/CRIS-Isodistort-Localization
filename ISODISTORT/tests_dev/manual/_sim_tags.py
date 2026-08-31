from isocore.backend.iso_wrapper import SubgroupInfo
from isocore.io.isodistort_cif import _irrep_kovalev_tag, _k_point_comment, _space_hm_symbol

# Simulate M3- at M for I4/mmm
sg = SubgroupInfo(
    index=0,
    space_group_number=129,
    space_group_symbol="P4/nmm",
    size=2,
    irrep_label="M3-",
    k_point_label="M",
    parent_sg=139,
    k_coordinates=["1", "1", "1"],
)
print("M3-", _irrep_kovalev_tag(sg), _k_point_comment(sg))

sg2 = SubgroupInfo(
    index=0,
    space_group_number=123,
    space_group_symbol="P4/mmm",
    size=12,
    irrep_label="LD1",
    k_point_label="LD",
    parent_sg=139,
    k_coordinates=["0", "0", "g"],
    k_parameters=["1/6"],
)
print("LD1", _irrep_kovalev_tag(sg2), _k_point_comment(sg2))

sg3 = SubgroupInfo(
    index=0,
    space_group_number=25,
    space_group_symbol="Pmm2",
    size=2,
    irrep_label="A1+",
    k_point_label="A",
    parent_sg=123,
    k_coordinates=["1/2", "1/2", "1/2"],
)
print("A1+ Nd", _irrep_kovalev_tag(sg3), _k_point_comment(sg3))
print("hm Pmm2", _space_hm_symbol("Pmm2"))
