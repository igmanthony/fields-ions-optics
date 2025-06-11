let array1 = Array(10000).fill(0)

let glob = {
    active_electrode: 1,
    emap_altered: false,
    dragging: false,
    dragging_frame: 0,
    initial_click: [-1, -1],
    brush: 1,
    electric_fields_calculated: false,
    voltages_were_updated: false,
    latest_ion_path: [],
    latest_ion: [],
    mouse_on_electrode_map: false,
    previous_mouse_electrode_position_x: -1,
    previous_mouse_electrode_position_y: -1,
    previous_mouse_electrode_position_x2: -1,
    previous_mouse_electrode_position_y2: -1,
    mouse_electrode_map_x: 0.0,
    mouse_electrode_map_y: 0.0,
    mouse_ef_map_x: 0.0,
    mouse_ef_map_y: 0.0,
    delta_ef_map_x: 0.0,
    delta_ef_map_y: 0.0,
    no_lines: true,
    splatlog: [],
    to_apply: [],
    rect: [0, 0],
    cellSize: 3,
    scale: 2,
    canvas_height: 600, // ususally 600
    canvas_width: 600, // ususally 600
    map_width: 200,
    map_height: 200,
    // undo
    history1: 0,
    history2: 0,
    history3: array1,
    // advanced randomization
    three_dimension_view: 0,
    efmap_altered: 0,
    gl: null,
    programInfo: null,
    bufferInfo: null,
    three_d_zoom_level: 0,
    efmap_pixels: null,
    fake_canvas: null,
    refresh_canvas: 1,
};

let state = {
    electrode_volts: [0, 500, 1000, -500, -1000],
    
    splattable: [-1, 0, 0, 0, 0], // -1 = background; 1 = splattable; 0 = not splattable
    dc_only: [1, 1, 1, 1, 1],
    magnet_mode: [0, 0, 0, 0, 0],
    
    pulse_starts: [0, 0, 0, 0, 0],
    pulse_ends: [0, 1000000, 1000000, 1000000, 1000000],
    
    electrode_frequencies: [0, 1.6e3, 1.6e3, 1.6e3, 1.6e3],
    electrode_offsets: [0, 0, 0, 0, 0],

    background_gas_mass: 28.0, // air
    pressure: 0.0, // in Pa

    current_mz: 100,
    current_x_position: 0.00,
    current_y_position: 0.00,
    current_x_velocity: 0.00,
    current_y_velocity: 0.00,
    random_ke: 0,

    max_time: 1000000, // steps
    fdm_accuracy: 10,
};

const turbo = [
    "#30123b", "#311542", "#32184a", "#341b51", "#351e58", "#36215f", "#372365", "#38266c",
    "#392972", "#3a2c79", "#3b2f7f", "#3c3285", "#3c358b", "#3d3791", "#3e3a96", "#3f3d9c",
    "#4040a1", "#4043a6", "#4145ab", "#4148b0", "#424bb5", "#434eba", "#4350be", "#4353c2",
    "#4456c7", "#4458cb", "#455bce", "#455ed2", "#4560d6", "#4563d9", "#4666dd", "#4668e0",
    "#466be3", "#466de6", "#4670e8", "#4673eb", "#4675ed", "#4678f0", "#467af2", "#467df4",
    "#467ff6", "#4682f8", "#4584f9", "#4587fb", "#4589fc", "#448cfd", "#438efd", "#4291fe",
    "#4193fe", "#4096fe", "#3f98fe", "#3e9bfe", "#3c9dfd", "#3ba0fc", "#39a2fc", "#38a5fb",
    "#36a8f9", "#34aaf8", "#33acf6", "#31aff5", "#2fb1f3", "#2db4f1", "#2bb6ef", "#2ab9ed",
    "#28bbeb", "#26bde9", "#25c0e6", "#23c2e4", "#21c4e1", "#20c6df", "#1ec9dc", "#1dcbda",
    "#1ccdd7", "#1bcfd4", "#1ad1d2", "#19d3cf", "#18d5cc", "#18d7ca", "#17d9c7", "#17dac4",
    "#17dcc2", "#17debf", "#18e0bd", "#18e1ba", "#19e3b8", "#1ae4b6", "#1be5b4", "#1de7b1",
    "#1ee8af", "#20e9ac", "#22eba9", "#24eca6", "#27eda3", "#29eea0", "#2cef9d", "#2ff09a",
    "#32f197", "#35f394", "#38f491", "#3bf48d", "#3ff58a", "#42f687", "#46f783", "#4af880",
    "#4df97c", "#51f979", "#55fa76", "#59fb72", "#5dfb6f", "#61fc6c", "#65fc68", "#69fd65",
    "#6dfd62", "#71fd5f", "#74fe5c", "#78fe59", "#7cfe56", "#80fe53", "#84fe50", "#87fe4d",
    "#8bfe4b", "#8efe48", "#92fe46", "#95fe44", "#98fe42", "#9bfd40", "#9efd3e", "#a1fc3d",
    "#a4fc3b", "#a6fb3a", "#a9fb39", "#acfa37", "#aef937", "#b1f836", "#b3f835", "#b6f735",
    "#b9f534", "#bbf434", "#bef334", "#c0f233", "#c3f133", "#c5ef33", "#c8ee33", "#caed33",
    "#cdeb34", "#cfea34", "#d1e834", "#d4e735", "#d6e535", "#d8e335", "#dae236", "#dde036",
    "#dfde36", "#e1dc37", "#e3da37", "#e5d838", "#e7d738", "#e8d538", "#ead339", "#ecd139",
    "#edcf39", "#efcd39", "#f0cb3a", "#f2c83a", "#f3c63a", "#f4c43a", "#f6c23a", "#f7c039",
    "#f8be39", "#f9bc39", "#f9ba38", "#fab737", "#fbb537", "#fbb336", "#fcb035", "#fcae34",
    "#fdab33", "#fda932", "#fda631", "#fda330", "#fea12f", "#fe9e2e", "#fe9b2d", "#fe982c",
    "#fd952b", "#fd9229", "#fd8f28", "#fd8c27", "#fc8926", "#fc8624", "#fb8323", "#fb8022",
    "#fa7d20", "#fa7a1f", "#f9771e", "#f8741c", "#f7711b", "#f76e1a", "#f66b18", "#f56817",
    "#f46516", "#f36315", "#f26014", "#f15d13", "#ef5a11", "#ee5810", "#ed550f", "#ec520e",
    "#ea500d", "#e94d0d", "#e84b0c", "#e6490b", "#e5460a", "#e3440a", "#e24209", "#e04008",
    "#de3e08", "#dd3c07", "#db3a07", "#d93806", "#d73606", "#d63405", "#d43205", "#d23005",
    "#d02f04", "#ce2d04", "#cb2b03", "#c92903", "#c72803", "#c52602", "#c32402", "#c02302",
    "#be2102", "#bb1f01", "#b91e01", "#b61c01", "#b41b01", "#b11901", "#ae1801", "#ac1601",
    "#a91501", "#a61401", "#a31201", "#a01101", "#9d1001", "#9a0e01", "#970d01", "#940c01",
    "#910b01", "#8e0a01", "#8b0901", "#870801", "#840701", "#810602", "#7d0502", "#7a0402",
];

const electrode_colors = [ "#ffffff", "#f44336", "#4CAF50", "#2196F3", "#795548", ];
const electrode_colors_100 = [ "#ffffff", "#FFCDD2", "#C8E6C9", "#BBDEFB", "#D7CCC8", ];
const electrode_colors_300 = [ "#ffffff", "#E57373", "#81C784", "#64B5F6", "#A1887F", ];

export {glob, state, turbo, electrode_colors, electrode_colors_100, electrode_colors_300}; 