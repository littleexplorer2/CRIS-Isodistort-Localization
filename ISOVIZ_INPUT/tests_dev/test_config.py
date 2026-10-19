from isoviz_input.backend.config_loader import CONFIG_PATH, get_config
from isoviz_input.backend.paths import DATA_DIR, INPUT_ROOT, ISOVIZ_ROOT, STRUCTURE_DIR


def test_settings_yaml_defines_input_paths() -> None:
    cfg = get_config()
    assert CONFIG_PATH.is_file()
    assert cfg.project_root == ISOVIZ_ROOT
    assert cfg.input_root == INPUT_ROOT
    assert cfg.data_dir == DATA_DIR
    assert cfg.structure_dir == STRUCTURE_DIR
    assert DATA_DIR.name == "data.csv"
    assert STRUCTURE_DIR.name == "subgroup.isoviz"
    assert cfg.best_model_folder_name == "Best_Model_Parameters"
    assert cfg.best_model_env_var == "BEST_MODEL_PARAMETERS_DIR"
    assert cfg.best_model_parent is None
    assert "ISOVIZ" in cfg.isoviz_env_vars
    assert "ISOViz.lnk" in cfg.launcher_names
    assert cfg.launcher_search_dirs == [
        cfg.project_root.parent / "launchers"
    ]
