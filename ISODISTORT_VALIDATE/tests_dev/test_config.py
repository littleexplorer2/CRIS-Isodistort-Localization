from isodistort_validate.config_loader import CONFIG_PATH, get_config
from isodistort_validate.compare_paths import COMPARE_ROOT, ITEM_DIR, TRUE_DIR, VALIDATE_ROOT


def test_settings_yaml_defines_compare_dirs() -> None:
    cfg = get_config()
    assert CONFIG_PATH.is_file()
    assert cfg.project_root == VALIDATE_ROOT
    assert cfg.item_dir == ITEM_DIR
    assert cfg.true_dir == TRUE_DIR
    assert ITEM_DIR.parent == COMPARE_ROOT
    assert ITEM_DIR.name == "item"
    assert TRUE_DIR.name == "true"
    assert cfg.pattern == "*.cif"
    assert cfg.lattice_tolerance == 1e-5
    assert cfg.coordinate_tolerance == 1e-5
    assert cfg.scalar_tolerance == 1e-5
    assert cfg.ignore_atom_order is False
    assert cfg.strict is False
