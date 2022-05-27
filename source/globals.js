import { turbo } from './constants.js';

let array1 = Array(10000).fill(0)

export default {
    electrode_volts: [0, 500, 1000, -500, -1000],
    splattable: [-1, 0, 0, 0, 0], // -1 = background; 1 = splattable; 0 = not splattable
    electrode_frequencies: [0, 1.6e3, 1.6e3, 1.6e3, 1.6e3],
    electrode_offsets: [0, 0, 0, 0, 0],
    pulse_starts: [0, 0, 0, 0, 0],
    pulse_ends: [0, 400000, 400000, 400000, 400000],
    active_electrode: 1,
    emap_altered: false,
    dragging: false,
    dragging_frame: 0,
    initial_click: [-1, -1],
    brush: 1,
    cmap: turbo,
    electric_fields_calculated: false,
    voltages_were_updated: false,
    latest_ion_path: [],
    latest_ion: [],
    
    current_mz: 100,
    current_x_position: 0.00,
    current_y_position: 0.00,
    current_x_velocity: 0.00,
    current_y_velocity: 0.00,

    mouse_on_electrode_map: false,
    previous_mouse_electrode_position_x: -1,
    previous_mouse_electrode_position_y: -1,
    previous_mouse_electrode_position_x2: -1,
    previous_mouse_electrode_position_y2: -1,
    mouse_electrode_map_x: 0.0,
    mouse_electrode_map_y: 0.0,
    no_lines: true,
    splatlog: [],
    dc_only: [1, 1, 1, 1, 1],
    magnet_mode: [0, 0, 0, 0, 0],

    to_apply: [],
    rect: [0, 0],
    pressure: 0.0, // in Pa
    background_gas_mass: 4.0, // helium

    cellSize: 3,
    scale: 2,
    map_width: 200,
    map_height: 200,
    max_time: 400000, // steps

    canvas_width: 600,
    canvas_height: 600,

    history1: 0,
    history2: 0,
    history3: array1,

};