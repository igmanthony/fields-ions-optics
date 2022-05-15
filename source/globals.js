import { turbo } from './constants.js';

export default {
    electrode_volts: [0, 500, 1000, -500, -1000],
    splattable: [-1, 0, 0, 0, 0], // -1 = background; 1 = splattable; 0 = not splattable
    electrode_frequencies: [0, 1.6e3, 1.6e3, 1.6e3, 1.6e3],
    electrode_offsets: [0, 0, 0, 0, 0],
    active_electrode: 1,
    emap_altered: false,
    dragging: false,
    dragging_frame: 0,
    initial_click: [-1, -1],
    brush: 1,
    cmap: turbo,
    electric_fields_calculated: false,
    latest_ion_path: [],
    latest_ion: [],
    
    current_mz: 100,
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
    dc_only: [true, true, true, true, true],

    // pressure
    pressure: 0.0,
    background_gas_mass: 4.0, // helium

};