import glob from './globals.js';
import { draw_electrodes, draw_over_electrodes, return_history, save_history } from './electrodes.js';
import { drawFieldLines, draw_electric_field, fly_ion } from './electric_fields.js';

import {
    electrode_colors,
    electrode_colors_300,
    default_splats,
    default_dcs,
    default_frequencies,
    default_pressure,
    default_offs,
    einzel_lens_values,
    quadrupole_values,
    quadrupole_frequencies,
    quadrupole_dcs,
    quadrupole_offsets,
    tof_values,
    tof_splats,
    funnel_values,
    funnel_frequencies,
    funnel_pressure,
    cyclotron_magnets,
    cyclotron_frequencies,
} from './constants.js';

function refreshToggles() {
    document.getElementById("splattable").disabled = !Boolean(glob.active_electrode);
    document.getElementById("splattable").checked = glob.splattable[glob.active_electrode];
    document.getElementById("dc_mode").disabled = !Boolean(glob.active_electrode);
    document.getElementById("dc_mode").checked = !glob.dc_only[glob.active_electrode];
    document.getElementById("magnet_mode").checked = glob.magnet_mode[glob.active_electrode];
    document.getElementById("sin_frequency").value = glob.electrode_frequencies[glob.active_electrode];
    document.getElementById("offset").value = glob.electrode_offsets[glob.active_electrode];
    if (document.getElementById("magnet_mode").checked) {
        document.getElementById("sin_frequency").disabled = true;
        document.getElementById("offset").disabled = true;
        document.getElementById("dc_mode").disabled = true;
        document.getElementById("splattable").disabled = true;
    } else if (!document.getElementById("dc_mode").checked) {
        document.getElementById("sin_frequency").disabled = true;
        document.getElementById("offset").disabled = true;
        document.getElementById("pulse_start").disabled = false;
        document.getElementById("pulse_end").disabled = false;
    } else {
        document.getElementById("sin_frequency").disabled = false;
        document.getElementById("offset").disabled = false;
        document.getElementById("pulse_start").disabled = true;
        document.getElementById("pulse_end").disabled = true;
    }
    document.getElementById("voltage_input").value = ((glob.active_electrode == 0) ? 0 : glob.electrode_volts[glob.active_electrode]);
}

function setupButtons(image, lens, tof, quadrupole, cyclotron, funnel) {
    const voltageInput = document.getElementById("voltage_input");
    voltageInput.oninput = () => {
        if (glob.active_electrode != 0) {
            const v = voltageInput.value;
            if (v <= 20000 && v >= -20000) {
                glob.electrode_volts[glob.active_electrode] = v;
                // voltageSlider.value = v;
                image.update_voltages(glob.electrode_volts);
                glob.voltages_were_updated = true;
            }
        }
    };

    // add box for first electrode on initial setup
    document.getElementById('electrode1').style.boxShadow = "0px 0px 3px 1px black";

    for (let i = 0; i < electrode_colors.length; i++) {
        document.getElementById(`electrode${i}`).style.background = electrode_colors_300[i];
        document.getElementById(`electrode${i}`).addEventListener("click", () => {
            glob.active_electrode = i;
            refreshToggles();
            for (let j = 0; j < electrode_colors.length; j++) { // electrode color buttons
                document.getElementById(`electrode${j}`)
                    .style.boxShadow = ((j == glob.active_electrode) ? "0px 0px 3px 1px black" : "");
            }
            voltageInput.value = ((i == 0) ? 0 : glob.electrode_volts[glob.active_electrode]);
        });
    }


    document.getElementById('brush1').style.boxShadow = "0px 0px 3px 1px black";
    for (let i = 1; i < 6; i++) {
        document.getElementById(`brush${i}`).addEventListener("click", (_) => {
            glob.brush = i;
            for (let i = 1; i < 6; i++) { // brush size buttons
                document.getElementById(`brush${i}`).style.boxShadow = ((i == glob.brush) ? "0px 0px 3px 1px black" : "");
            }
        });
    }

    document.getElementById('brush6').addEventListener("click", (event) => {
        return_history(image);
    });

    document.getElementById("splattable").addEventListener("change", (event) => { 
        glob.splattable[glob.active_electrode] = (event.target.checked && glob.active_electrode != 0)? 1 : 0; // can't pass bools so this suffices
        draw_electrodes(image);
    });

    document.getElementById("dc_mode").addEventListener("change", (event) => { 
        if (event.target.checked && glob.active_electrode != 0) {
            glob.dc_only[glob.active_electrode] = 0; // can't pass bools so this suffices
            document.getElementById("sin_frequency").disabled = false;
            document.getElementById("offset").disabled = false;
            document.getElementById("pulse_start").disabled = true;
            document.getElementById("pulse_end").disabled = true;
        } else {
            glob.dc_only[glob.active_electrode] = 1; // can't pass bools so this suffices
            document.getElementById("sin_frequency").disabled = true;
            document.getElementById("offset").disabled = true;
            document.getElementById("pulse_start").disabled = false;
            document.getElementById("pulse_end").disabled = false;
        }
    });

    document.getElementById("magnet_mode").addEventListener("change", (event) => { 
        if (event.target.checked && glob.active_electrode != 0) {
            glob.magnet_mode[glob.active_electrode] = 1;
            document.getElementById("sin_frequency").disabled = true;
            document.getElementById("offset").disabled = true;
            document.getElementById("dc_mode").disabled = true;
            document.getElementById("splattable").disabled = true;
            document.getElementById("pulse_start").disabled = true;
            document.getElementById("pulse_end").disabled = true;
        } else {
            glob.magnet_mode[glob.active_electrode] = 0;
            document.getElementById("sin_frequency").disabled = false;
            document.getElementById("offset").disabled = false;
            document.getElementById("dc_mode").disabled = false;
            document.getElementById("splattable").disabled = false;
            document.getElementById("pulse_start").disabled = false;
            document.getElementById("pulse_end").disabled = false;
        }
        refreshToggles();
        draw_electrodes(image);
    });

    document.getElementById("update").addEventListener("click", (_) => {
        image.update_voltages(glob.electrode_volts);
        image.update_solids(glob.splattable);
        image.update_magnets(glob.magnet_mode);
        refreshToggles();
        console.log(`Updated voltages ${glob.electrode_volts} and solids ${glob.splattable}`);
        image.update_timeline(
            glob.dc_only,
            glob.electrode_frequencies,
            glob.pulse_starts,
            glob.pulse_ends
        );
        console.log("Made timeline");
        image.generate_electric_fields();
        console.log("Updated electrodes & generated electric fields");
        glob.electric_fields_calculated = true;
        draw_electric_field(image);
    })

    document.getElementById("fly'm").addEventListener("click", (_) => {
        glob.current_x_position = parseFloat(document.getElementById("x-position").value);
        glob.current_y_position = glob.map_height - parseFloat(document.getElementById("y-position").value);
        console.log(glob.current_x_position);
        console.log(glob.current_y_position);
        fly_ion(image, glob.current_x_position, glob.current_y_position);
    })

    // TOP
    document.getElementById("clear").addEventListener("click", (_) => {
        glob.no_lines = true;
        image.clear();
        refreshToggles();
        save_history(image);
        draw_electrodes(image);
    })

    document.getElementById("einzel_lens").addEventListener("click", (_) => setElectrodeValues(image, lens, einzel_lens_values, default_splats, default_frequencies, default_dcs, default_pressure, default_offs, default_offs));
    document.getElementById("tof").addEventListener("click", (_) => setElectrodeValues(image, tof, tof_values, tof_splats, default_frequencies, default_dcs, default_pressure, default_offs, default_offs));
    document.getElementById("quadrupole").addEventListener("click", (_) => setElectrodeValues(image, quadrupole, quadrupole_values, default_splats, quadrupole_frequencies, quadrupole_dcs, default_pressure, quadrupole_offsets, default_offs));
    document.getElementById("cyclotron").addEventListener("click", (_) => setElectrodeValues(image, cyclotron, quadrupole_values, default_splats, cyclotron_frequencies, quadrupole_dcs, default_pressure, default_offs, cyclotron_magnets));
    document.getElementById("ion_funnel").addEventListener("click", (_) => setElectrodeValues(image, funnel, funnel_values, default_splats, funnel_frequencies, quadrupole_dcs, funnel_pressure, default_offs, default_offs));

    document.getElementById("SIMION_SAVE").addEventListener("click", (_) => {
        image.generate_electrodes();
        let blob = image.save_simion_pa();
        let filename = "electrodes.pa";
        var file = new Blob([blob], { type: 'octet/stream' });
        if (window.navigator.msSaveOrOpenBlob) {// IE10+
            window.navigator.msSaveOrOpenBlob(file, filename);
        } else { // Others
            let a = document.createElement("a"),
                url = URL.createObjectURL(file);
            a.href = url;
            a.download = filename;
            document.body.appendChild(a);
            a.click();
            setTimeout(
                function () {
                    document.body.removeChild(a);
                    window.URL.revokeObjectURL(url);
                },
                0
            );
        }
    });

    // advanced
    document.getElementById("field_lines").addEventListener("click", (_) => { drawFieldLines(image); }); // can I remove the closure?
    document.getElementById("mz_ratio").addEventListener("input", (_) => { glob.current_mz = parseFloat(document.getElementById("mz_ratio").value); });
    document.getElementById("ion_x_velocity").addEventListener("input", (_) => { glob.current_x_velocity = parseFloat(document.getElementById("ion_x_velocity").value); });
    document.getElementById("ion_y_velocity").addEventListener("input", (_) => { glob.current_y_velocity = parseFloat(document.getElementById("ion_y_velocity").value); });
    
    
    document.getElementById("fdm_threshold").addEventListener("change", (_) => {
        image.update_fdm_threshold(parseFloat(document.getElementById("fdm_threshold").value));
    });
    
    
    document.getElementById("pressure").addEventListener("input", (_) => {
        glob.pressure = parseFloat(document.getElementById("pressure").value);
        image.update_pressure(glob.pressure);
    });
    document.getElementById("bg_gas").addEventListener("change", (_) => {
        glob.background_gas_mass = parseFloat(document.getElementById("bg_gas").value);
        console.log(glob.background_gas_mass);
        image.update_gas_mass(glob.background_gas_mass);
    });
    document.getElementById("sin_frequency").addEventListener("input", (a) => {
        glob.electrode_frequencies[glob.active_electrode] = parseFloat(document.getElementById("sin_frequency").value);
    });
    document.getElementById("offset").addEventListener("input", (a) => {
        glob.electrode_offsets[glob.active_electrode] = parseFloat(document.getElementById("offset").value);
        image.update_offsets(glob.electrode_offsets);
    });
    document.getElementById("pulse_start").addEventListener("input", (_) => {
        glob.pulse_starts[glob.active_electrode] = parseFloat(document.getElementById("pulse_start").value);
    });
    document.getElementById("pulse_end").addEventListener("input", (_) => {
        glob.pulse_ends[glob.active_electrode] = parseFloat(document.getElementById("pulse_end").value);
    });
    
}

function setElectrodeValues(image, map, volts, splattable, frequencies, dc_only_settings, pressure, offsets, magnets) {
    let width = glob.map_width / 2;
    let height = glob.map_height / 2;

    for (let x=0; x < width; x++)  {
        for (let y=0; y < height; y++) {
            let color = map[(x + width * y)];
            image.brush(x, y, color);
        }
    }
    glob.emap_altered = true;
    for (let i=0; i < volts.length; i++) {
        glob.electrode_volts[i] = volts[i]; // no idea why this works to reset the values and just setting electrode_volts equal to lens values doesn't work
        glob.splattable[i] = splattable[i]; 
    }
    for (let i = 0; i < frequencies.length; i++) {
        glob.electrode_frequencies[i] = frequencies[i];
        glob.dc_only[i] = dc_only_settings[i];
        glob.electrode_offsets[i] = offsets[i];
        glob.magnet_mode[i] = magnets[i];
    }
    image.update_pressure(pressure);
    image.update_gas_mass(28.0);
    image.update_offsets(offsets);
    image.update_magnets(glob.magnet_mode);
    document.getElementById("pressure").value = pressure;
    document.getElementById('bg_gas').getElementsByTagName('option')[0].selected = 'selected'
    refreshToggles();
    save_history(image)
}

export { setupButtons };