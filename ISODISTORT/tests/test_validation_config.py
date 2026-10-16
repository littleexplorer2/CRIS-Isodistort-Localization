from tests.validation.backend.config import CONFIG_PATH, get_config
from tests.validation.backend.paths.service import COMPARE_ROOT, ITEM_DIR, TRUE_DIR, VALIDATE_ROOT
from tests.validation.backend.config import service as config_loader


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


def test_config_path_override_keeps_relative_paths_beside_selected_file(tmp_path, monkeypatch) -> None:
    selected = tmp_path / "custom-config" / "settings.yaml"
    monkeypatch.setattr(config_loader, "CONFIG_PATH", selected)

    assert get_config().resolve_path("../custom-compare") == tmp_path / "custom-compare"
