"""
English UI strings for the web page and API console output.

Keys use scene.action names such as load.done. Placeholders use Python
str.format syntax ({n}, {sg}) or numbered web placeholders ({0}).
"""

MESSAGES: dict[str, str] = {
    'load.done': '[Loaded] Space group #{sg} ({sym}), {n} atoms',
    'subgroups.found': '[Subgroup enumeration] {n} isotropy subgroups found',
    'subgroups.more': '  ... {n} more',
    'subgroup.not_found': 'Subgroup index {idx} not found',
    'path.selected': '[Path selected] {desc}',
    'modes.found': '[Mode calculation] {n} distortion modes',
    'mode.sites': '  {irrep:<8s} {opd:<6s} affecting {n} Wyckoff position(s)',
    'mode.empty': (
        "This subgroup has no displacement modes on the structure's Wyckoff positions "
        "(type filtering is applied at the mode-calculation stage)"
    ),
    'distortion.generated': '[Distortion] mode {irrep}, amplitude {amp}, atoms {n1} -> {n2} (supercell factor {r:g})',
    'mode.invalid': 'Mode {label} does not exist',
    'export.default': '[Auto-export] CIF written: {path}',
    'export.done': '[Export done] {n} file(s):',
    'domains.found': '[Domains] {n} domain(s) generated',
    'method1.result': '[Method 1] {n} candidate(s) after filtering',
    'method2.result': '[Method 2] Subgroup #{idx}: {n} mode(s)',
    'method3.result': '[Method 3] {n} candidate(s) from constrained search',
    'method4.result': '[Method 4] Decomposition done: {n} mode(s), RMS residual = {rms:.6f}',
    'err.load_first': 'Load a structure first (load_structure)',
    'err.select_path_first': 'Select a phase path first (select_path) or run Method 2',
    'err.generate_first': 'Generate a distorted structure first (generate_distortion)',
    'err.mode_requires_path': 'Compute modes first via select_path or search_method_2',
    'err.domains_need_path': 'Select a phase path first (select_path)',
    'err.domains_not_in_list': "The path's subgroup is not in the candidate list",
    'hStatus': 'Session',
    'hCif': 'Parent CIF',
    'hTypes': 'Distortion Types',
    'hDist': 'Distortion Page',
    'btn.load': 'Load',
    'btn.apply': 'Apply',
    'btn.m1': 'Run Method 1',
    'btn.m2': 'Run Method 2',
    'btn.m3': 'Run Method 3',
    'btn.m4': 'Run Method 4',
    'btn.subs': 'List subgroups',
    'btn.exp': 'Export',
    'l.sg': 'SG number',
    'l.max': 'Maximal subgroups only',
    'l.idx': 'Subgroup idx',
    's.direct': 'Direct k-point search',
    'l.gen': 'Generate DB if missing',
    'l.sg3': 'Space group type',
    'l.basis': 'Choose a representative basis:',
    'l.match': 'Matching',
    'l.thresh': 'Threshold',
    'l.fmt': 'Formats',
    'st.noStruct': 'No parent structure loaded',
    'st.loaded': 'SG #{0} ({1}), {2} atoms',
    'st.types': 'Types',
    'st.subgroups': 'Subgroups',
    'st.modes': 'Modes',
    'st.distorted': 'Distorted atoms',
    'st.wait': 'Computing, please wait…',
    'st.elapsed': 'Elapsed {0}',
    'st.busyHint': 'Local engine is running (WSL / iso). Keep this tab open — the page will update when finished.',
    'st.busyStatus': 'Busy: {0} · {1}',
    'st.noModes': 'No modes; run Method 2 first',
    'ok.loaded': 'Parent structure loaded',
    'ok.generated': 'Distortion generated, {0} atoms',
    'err.noFile': 'Choose a file first',
    'err.needIdx': 'Enter a subgroup index',
    'err.noModes': 'No modes available',
    'err.noExport': 'Nothing to export yet',
    'hIdx': 'idx',
    'hK': 'k point',
    'hIrrep': 'Irrep',
    'hRouteStatus': 'route status',
    'hKnownRoutes': 'known routes',
    'hDir': 'Dir',
    'hKActive': 'k-active',
    'hNrep': 'Reps',
    'hMode': 'Mode',
    'hAs': 'As (angstrom)',
    'hAp': 'Ap (angstrom)',
    'hRawCoefficient': 'raw coefficient',
    'hNormfactor': 'normfactor (1/angstrom)',
    'hGen': 'Generator',
    'web.title': 'ISODISTORT: search',
    'types.title': 'Types of distortions to be considered',
    'btn.ok': 'OK',
    'btn.change': 'Change',
    'btn.stop': 'Stop',
    'type.note': (
        'Important: You must click on Change to implement any changes in the '
        'above type of distortions to be considered.'
    ),
    'ok.types': 'Types updated.',
    'l.cs': 'Crystal system(s):',
    'l.sgsel': 'Space-group symmetry:',
    'l.convLat': 'Conventional lattice:',
    'l.primLat': 'Primitive lattice:',
    'm1.sgHint': '({0} reachable subgroup space groups, filtered by the parent structure)',
    'm1.latNote': '',
    'err.badNsup': 'Superposed IR count must be a positive integer.',
    'ok.nsup': 'Superposed IR count set to {0}.',
    'l.noChoice': 'no choice',
    'l.sublat': 'Direct sublattice (a,b,c) (default 1,1,1):',
    'm1.latHint': '(e.g. 2,2,2 = double each axis)',
    'm1.title': 'Method 1: Search over all special k points',
    'm2.title': 'Method 2: General method - search over specific k points',
    'm3.title': 'Method 3: Search over arbitrary k points for specified space group and lattice',
    'm4.title': 'Method 4: Mode decomposition of a distorted structure',
    'l.kp': 'Specify k point:',
    'l.kvec': '<i>k</i> vector {0}:',
    'l.ir': 'Irreducible representation (IR):',
    'l.nmod': '# of independent incommensurate modulations:',
    'nmod.hint': '(0 = all folding k; n≥1 = harmonics of this Method 2 q, plus Gamma)',
    'm2.nmodNote': (
        'Important: You must click on Change to implement any changes in the number of '
        'independent incommensurate modulations. 0 = 3D lock-in (every parent k that '
        'folds into the child cell, including secondary IRs such as GM/M). '
        'n≥1 keeps harmonics of the single q selected for this subgroup, plus Gamma. '
        'Method 2 has one primary q, so 1, 2 and 3 do not add a second modulation.'
    ),
    'err.badNmod': 'Number of independent modulations must be an integer 0–3.',
    'ok.nmod': 'Independent modulations set to {0}.',
    'l.sgsel3': 'Select either space group symmetry:',
    'l.pg': 'or point group (crystal class):',
    'l.lattice': 'Specify a real-space sublattice of the parent lattice with',
    'l.cent': 'centering:',
    'l.m4file': 'Upload distorted structure from CIF file:',
    'l.m4matching': 'Atom-matching method:',
    'l.m4threshold': 'Robust distance threshold (angstrom):',
    'l.m4origin': 'Known origin shift (daughter fractional x,y,z; optional):',
    'm4.nearest': 'nearest-site (global minimum assignment)',
    'm4.robust': 'robust (reject assignments beyond threshold)',
    'm4.residualSummary': 'RMS residual {0} angstrom · max |residual component| {1} angstrom',
    'm4.strainSummary': 'Homogeneous strain ({0}): {1}',
    'st.stopping': 'Stopping the server and releasing the port…',
    'st.srvDown': 'Cannot reach the server — it may have stopped. Close this page or rerun main_web.py.',
    'st.noSubs': 'No subgroups.',
    'err.noCif': 'Load a parent CIF first.',
    'err.badSublat': 'Invalid direct sublattice; expected a,b,c of positive integers.',
    'ok.subs': '{0} subgroup(s) found.',
    'm1.orderParam': 'Order parameter:',
    'web.reading': 'Reading CIF file...',
    'web.done': 'Done.',
    'web.lattice': 'Lattice parameters:',
    'web.prefs': 'Default space-group preferences:',
    'dist.note': (
        'Methods 1–4 only compute results. This panel downloads them: '
        'filtered tables (all Methods) and subgroup structure files '
        '(Methods 1–3 ZIP).'
    ),
    'scope.all': 'all',
    'scope.none': 'none',
    'st.scope': 'Scope',
    'l.opd': 'Order parameter direction (OPD):',
    'l.nsup': 'Change number of superposed IRs:',
    'm2.nsup_note': 'Important: You must click on Change to implement any changes in the number of superposed IRs.',
    'm2.paramKpSelected': (
        'Note: you selected a parametric (incommensurate) k point. '
        'The local engine enumerates isotropy subgroups and computes complete '
        'displacive modes via smodes + the child space group (3D lock-in / (3+d) '
        'harmonics, controlled by nmod).'
    ),
    'm2.enumKp': 'Enumerating subgroups over all irreps of k point {0} (k group {1}/{2})...',
    'm2.noSubsAtKp': (
        'The local engine could not enumerate/generate subgroups for this k point '
        '(common for parametric k points such as LD/DT).'
    ),
    'm2.chooseNext': 'Choose next step:',
    'm2.localCompute': '① Compute with local resources',
    'm2.localComputeDesc': (
        "Use your machine to generate this k point's subgroup database "
        '(Generate isotropy subgroups; may take minutes to hours, then cached)'
    ),
    'm2.gotoOfficial': '② Retry on the ISODISTORT website',
    'm2.gotoOfficialDesc': (
        'Run Method 2 with the same parent CIF and (k point, parameters) on the '
        'website (this k point generates subgroups there)'
    ),
    'm2.subsFound': (
        'Enumerated {0} subgroup(s). Click a row to view its mode basis '
        '(official order parameter direction page).'
    ),
    'm3.subsFound': (
        'Found {0} local Method 3 embedding candidate(s). '
        'Click a row to view the representative known route.'
    ),
    'm3.emptyHint': (
        'No isotropy subgroups matched these Method 3 constraints. '
        'Check space-group / point-group filters and the supercell basis. '
        'For commensurate supercells (e.g. 1×1×6 → LD g=1/6), enable '
        '“Generate isotropy subgroups database if missing” under Method 2 '
        'if the parametric-k database is not cached yet. '
        'Default and P/A/B/C/I/F/R centerings are interpreted as exact primitive '
        'translation lattices. Reciprocal search is not supported locally.'
    ),
    'm2.filter': 'Filter:',
    'm2.clearFilter': 'Clear',
    'm2.showFiltered': 'Show filtered rows only',
    'm2.downloadTxt': 'Download filtered (txt)',
    'm2.downloadCsv': 'Download filtered (csv)',
    'm2.filteredCount': '{0} / {1} after filtering',
    'm2.paramKNote': (
        'This subgroup belongs to a parametric k point. Complete displacive modes '
        'are computed with smodes and the child space-group identity representation '
        '(nmod=0 lock-in includes folding harmonics and secondary IRs).'
    ),
    'lGenDb': 'Generate isotropy subgroups database if missing',
    'm2.genDbHelp': (
        'If this parametric k point has no prebuilt isotropy-subgroup list in the local '
        'ISOTROPY data files, enabling generation lets your machine run iso to create and '
        'cache that list (often minutes to hours the first time; later searches reuse the '
        'cache until you delete selected entries via Manage cached subgroup databases).'
    ),
    'm2.genDbWarn': (
        'Generation can take minutes to hours; cached files are stored in the local WSL '
        'staging temp directory as i*.iso (not resources/isobyu/) and can be listed or deleted in '
        'Manage cached subgroup databases.'
    ),
    'm2.genDbManage': 'Manage cached subgroup databases',
    'm2.genDbManageHide': 'Hide cache list',
    'm2.genDbEmpty': 'No generated subgroup databases are cached yet.',
    'm2.genDbSelectAll': 'Select all',
    'm2.genDbSelectNone': 'Select none',
    'm2.genDbDelete': 'Delete selected',
    'm2.genDbDeleted': 'Deleted {0} file(s).',
    'm2.genDbColName': 'File',
    'm2.genDbColSg': 'Parent SG',
    'm2.genDbColIr': 'Irrep',
    'm2.genDbColK': 'kparam',
    'm2.genDbColSize': 'Size',
    'm2.genDbColTime': 'Saved',
    'm2.genDb': '(generating subgroup database)',
    'm2.genDbRetry': 'Retry with local database generation',
    'm3.direct': 'Specify a real-space sublattice of the parent lattice with',
    'm3.reciprocal': 'Specify a primitive reciprocal-space superlattice',
    'm3.centDefault': 'Default',
    'err.reciprocal': 'The local engine does not support reciprocal (reciprocal-space supercell) mode; use direct.',
    'err.occSupercell': 'Occupational modes require the same subgroup supercell basis used at generation.',
    'occ.note': '(approximate: failed the subgroup symmetry check)',
    'dist.title': 'Distortion',
    'dist.downloadAll': 'Download all files (ZIP)',
    'dist.fmtLabel': 'Download formats:',
    'dist.fmtCif': 'CIF file',
    'dist.fmtIsoviz': 'Save interactive distortion',
    'dist.fmtModes': 'Complete modes details',
    'dist.fmtTopas': 'TOPAS.STR',
    'dist.noFmt': 'Select at least one export format.',
    'dist.methodLabel': 'Export source (exactly one Method):',
    'dist.method1': 'Method 1',
    'dist.method2': 'Method 2',
    'dist.method3': 'Method 3',
    'dist.method4': 'Method 4',
    'dist.tableLabel': 'Filtered result table:',
    'dist.tableNote': (
        'Table and ZIP downloads both use that Method’s current filters '
        '(all matching rows, even if the table is set to show filtered rows '
        'only). ZIP contains only the selected file formats. Method 4 has a '
        'table only — no ZIP.'
    ),
    'dist.noTable': 'The selected Method has no result table yet; run that Method first.',
    'dist.zipNotM4': (
        'Method 4 is a mode-amplitude table, not a subgroup list. Use '
        'Download filtered (txt/csv). ZIP is only for Methods 1–3.'
    ),
    'dist.noMethod': 'The selected Method has no subgroups to export; run that Method first.',
    'dist.zipWait': (
        'Building ZIP… keep this tab open. Non-CIF formats re-run Method 2 per '
        'subgroup to fill modes (can take a while). Parametric k points (LD/DT, …) '
        'use smodes/(3+d) complete modes (nmod). When strain is selected, each '
        'subgroup also resolves and validates its canonical ISO rank-[12] strain '
        'basis. Some special-k paths may omit Wyckoff sites or return an empty '
        'BUSH table versus the website.'
    ),
    'dist.zipParamNote': (
        'Note: one or more selected subgroups use a parametric k point. '
        'Displacement modes are filled with the local smodes/(3+d) complete-mode '
        'engine (nmod); selected strain modes use the same validated canonical '
        'rank-[12] export path as other subgroups.'
    ),
    'dist.paramFormatDefault': (
        'This Method includes a parametric k point. ZIP formats other than CIF '
        'now include smodes/(3+d) complete displacive modes (nmod=0 lock-in by '
        'default). Selected strain modes are included after canonical fixed-space '
        'validation.'
    ),
    'dist.zipDone': 'Downloaded ZIP ({0} subgroup(s)).',
    'dist.zipFail': 'ZIP download failed.',
    'hPrefs': 'Space-Group Preferences',
    'prefs.note': (
        'The local engine uses these fixed defaults (international standard setting, '
        'i.e. the website defaults); they cannot be modified locally.'
    ),
    'prefs.monoAxes': 'Monoclinic axes:',
    'prefs.monoCell': 'Monoclinic cell choice:',
    'prefs.orthoAxes': 'Orthorhombic axes:',
    'prefs.trigAxes': 'Trigonal axes:',
    'prefs.origin': 'Origin choice:',
    'prefs.ssg': 'Superspace group setting:',
    'prefs.ssgValue': 'standard (IT-C)',
    'prefs.fixedNote': (
        'Reason: the local iso binary only supports the international standard '
        'setting; custom settings (axes/cell choice/origin/SSG etc.) cause a Syntax '
        "error, so the website's preference-editing panel is not available locally; "
        'computations always use the defaults above. Superspace-group setting is '
        'fixed to standard (IT-C).'
    ),
    'm3.centEnd': 'centering.',
    'm1.wait': "Enumerating subgroups over all special k points (like the website's database query); please wait...",
}
