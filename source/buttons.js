import glob from './globals.js';
import { draw_electrodes, draw_over_electrodes, makeTimeline } from './electrodes.js';
import { drawFieldLines, draw_electric_field, fly_ion } from './electric_fields.js';

import {
    electrode_colors,
    electrode_colors_300,
    default_splats,
    default_dcs,
    default_frequencies,
    default_pressure,
    default_offsets,
    einzel_lens,
    einzel_lens_values,
    quadrupole,
    quadrupole_values,
    quadrupole_frequencies,
    quadrupole_dcs,
    quadrupole_offsets,
    tof,
    tof_values,
    tof_splats,
    funnel,
    funnel_values,
    funnel_frequencies,
    funnel_pressure,
} from './constants.js';

function refreshToggles() {
    document.getElementById("splattable").disabled = !Boolean(glob.active_electrode);
    document.getElementById("splattable").checked = glob.splattable[glob.active_electrode];
    document.getElementById("dc_mode").disabled = !Boolean(glob.active_electrode);
    document.getElementById("dc_mode").checked = !glob.dc_only[glob.active_electrode];
    document.getElementById("sin_frequency").value = glob.electrode_frequencies[glob.active_electrode];
    document.getElementById("offset").value = glob.electrode_offsets[glob.active_electrode];
    if (!document.getElementById("dc_mode").checked) {
        console.log("disabling freqoffset");
        document.getElementById("sin_frequency").disabled = true;
        document.getElementById("offset").disabled = true;
    } else {
        console.log("enabling freqoffset");
        document.getElementById("sin_frequency").disabled = false;
        document.getElementById("offset").disabled = false;
    }
}

function setupButtons(image) {
    // voltage slider and input
    // const voltageSlider = document.getElementById("voltage_slider");
    const voltageInput = document.getElementById("voltage_input");
    // voltageSlider.oninput = () => {
    //     if (glob.active_electrode != 0) {
    //         glob.electrode_volts[glob.active_electrode] = voltageSlider.value;
    //         voltageInput.value = voltageSlider.value;
    //     }
    // };
    voltageInput.oninput = () => {
        if (glob.active_electrode != 0) {
            const v = voltageInput.value;
            if (v <= 20000 && v >= -20000) {
                glob.electrode_volts[glob.active_electrode] = v;
                // voltageSlider.value = v;
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
            // voltageSlider.value = ((i == 0) ? 0 : glob.electrode_volts[glob.active_electrode]);
            voltageInput.value = ((i == 0) ? 0 : glob.electrode_volts[glob.active_electrode]);
        });
    }


    document.getElementById('brush1').style.boxShadow = "0px 0px 3px 1px black";
    for (let i = 1; i < 5; i++) {
        document.getElementById(`brush${i}`).addEventListener("click", (_) => {
            glob.brush = i;
            for (let i = 1; i < 5; i++) { // brush size buttons
                document.getElementById(`brush${i}`).style.boxShadow = ((i == glob.brush) ? "0px 0px 3px 1px black" : "");
            }
        });
    }
    document.getElementById("splattable").addEventListener("change", (event) => { 
        if (event.target.checked && glob.active_electrode != 0) {
            glob.splattable[glob.active_electrode] = 1; // can't pass bools so this suffices
        } else {
            glob.splattable[glob.active_electrode] = 0; // can't pass bools so this suffices
        }
    });

    document.getElementById("dc_mode").addEventListener("change", (event) => { 
        if (event.target.checked && glob.active_electrode != 0) {
            glob.dc_only[glob.active_electrode] = 0; // can't pass bools so this suffices
            document.getElementById("sin_frequency").disabled = false;
            document.getElementById("offset").disabled = false;
        } else {
            glob.dc_only[glob.active_electrode] = 1; // can't pass bools so this suffices
            document.getElementById("sin_frequency").disabled = true;
            document.getElementById("offset").disabled = true;
        }
    });


    document.getElementById("update").addEventListener("click", (_) => {
        image.update_voltages(glob.electrode_volts);
        image.update_splattables(glob.splattable);
        refreshToggles();
        console.log(`Updated voltages ${glob.electrode_volts} and splattables ${glob.splattable}`);
        image.update_timeline(makeTimeline(image.timeline_length()));
        console.log("Made timeline");
        image.generate_electric_fields();
        console.log("Updated electrodes & generated electric fields");
        glob.electric_fields_calculated = true;
        draw_electric_field(image);
    })

    document.getElementById("fly'm").addEventListener("click", (_) => {
        fly_ion(image, glob.current_x_position, glob.current_y_position);
    })

    // TOP
    document.getElementById("clear").addEventListener("click", (_) => {
        glob.no_lines = true;
        image.clear();
        refreshToggles();
        draw_electrodes(image);
    })

 


    document.getElementById("einzel_lens").addEventListener("click", (_) => setElectrodeValues(image, einzel_lens, einzel_lens_values, default_splats, default_frequencies, default_dcs, default_pressure, default_offsets));
    document.getElementById("tof").addEventListener("click", (_) => setElectrodeValues(image, tof, tof_values, tof_splats, default_frequencies, default_dcs, default_pressure, default_offsets));
    document.getElementById("quadrupole").addEventListener("click", (_) => setElectrodeValues(image, quadrupole, quadrupole_values, default_splats, quadrupole_frequencies, quadrupole_dcs, default_pressure, quadrupole_offsets));
    // document.getElementById("cyclotron").addEventListener("click", (_) => setElectrodeValues(quadrupole, quadrupole_values, default_splats, quadrupole_frequencies, quadrupole_dcs));
    document.getElementById("ion_funnel").addEventListener("click", (_) => setElectrodeValues(image, funnel, funnel_values, default_splats, funnel_frequencies, quadrupole_dcs, funnel_pressure, default_offsets));

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

function setElectrodeValues(image, map, volts, splattable, frequencies, dc_only_settings, pressure, offsets) {
    let dimension = Math.sqrt(map.length); 
    for (let x=0; x < dimension; x++)  {
        for (let y=0; y < dimension; y++) {
            let color = map[(x + dimension * y)];
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
    }
    image.update_pressure(pressure);
    image.update_gas_mass(28.0);
    image.update_offsets(offsets);
    document.getElementById("pressure").value = pressure;
    document.getElementById('bg_gas').getElementsByTagName('option')[0].selected = 'selected'
    refreshToggles();
}

export { setupButtons };