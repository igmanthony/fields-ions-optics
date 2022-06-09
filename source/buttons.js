import { glob, state, electrode_colors, electrode_colors_300 } from './globals.js';
import { draw_electrodes, save_history } from './electrodes.js';
import { draw_electric_field, fly_ion, setupElectricFieldMap } from './fields.js';

function setupButtons(image) { //, lens, tof, quadrupole, cyclotron, funnel) {
    const voltageInput = document.getElementById("voltage_input");
    voltageInput.oninput = () => {
        if (glob.active_electrode != 0) {
            const v = voltageInput.value;
            if (v <= 20000 && v >= -20000) {
                state.electrode_volts[glob.active_electrode] = v;
                image.update_voltages(state.electrode_volts);
                glob.voltages_were_updated = true;
            }
            glob.efmap_altered = 1;
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
            voltageInput.value = ((i == 0) ? 0 : state.electrode_volts[glob.active_electrode]);
        });
    }


    document.getElementById('brush1').style.boxShadow = "0px 0px 3px 1px black";
    for (let i = 1; i < 6; i++) {
        document.getElementById(`brush${i}`).addEventListener("click", (_) => {
            glob.brush = i;
            for (let j = 1; j < 6; j++) { // brush size buttons
                document.getElementById(`brush${j}`).style.boxShadow = ((j == glob.brush) ? "0px 0px 3px 1px black" : "");
            }
        });
    }
    document.getElementById('brush6').addEventListener("click", (_) => { return_history(image); });

    document.getElementById("splattable").addEventListener("change", (event) => { 
        state.splattable[glob.active_electrode] = (event.target.checked && glob.active_electrode != 0)? 1 : 0; // can't pass bools so this suffices
        draw_electrodes(image);
    });

    document.getElementById("dc_mode").addEventListener("change", (event) => { 
        if (event.target.checked && glob.active_electrode != 0) {
            state.dc_only[glob.active_electrode] = 0; // can't pass bools so this suffices
            document.getElementById("sin_frequency").disabled = false;
            document.getElementById("offset").disabled = false;
            document.getElementById("pulse_start").disabled = true;
            document.getElementById("pulse_end").disabled = true;
        } else {
            state.dc_only[glob.active_electrode] = 1; // can't pass bools so this suffices
            document.getElementById("sin_frequency").disabled = true;
            document.getElementById("offset").disabled = true;
            document.getElementById("pulse_start").disabled = false;
            document.getElementById("pulse_end").disabled = false;
        }
    });

    document.getElementById("magnet_mode").addEventListener("change", (event) => { 
        if (event.target.checked && glob.active_electrode != 0) {
            state.magnet_mode[glob.active_electrode] = 1;
            document.getElementById("sin_frequency").disabled = true;
            document.getElementById("offset").disabled = true;
            document.getElementById("dc_mode").disabled = true;
            document.getElementById("splattable").disabled = true;
            document.getElementById("pulse_start").disabled = true;
            document.getElementById("pulse_end").disabled = true;
        } else {
            state.magnet_mode[glob.active_electrode] = 0;
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
        update_stuff(image);
        console.log(`Updated voltages ${state.electrode_volts} and solids ${state.splattable}`);
        image.update_timeline(
            state.dc_only,
            state.electrode_frequencies,
            state.pulse_starts,
            state.pulse_ends
        );
        console.log("Made timeline");
        image.generate_electric_fields();
        glob.electric_fields_calculated = true;
        glob.fake_canvas = null;
        draw_electric_field(image);
    });


    document.getElementById("fly'm").addEventListener("click", (_) => {
        state.current_x_position = parseFloat(document.getElementById("x-position").value);
        state.current_y_position = glob.map_height - parseFloat(document.getElementById("y-position").value);
        fly_ion(image, state.current_x_position, state.current_y_position);
    });

    // TOP
    document.getElementById("clear").addEventListener("click", (_) => {
        glob.no_lines = true;
        image.clear();
        refreshToggles();
        save_history(image);
        draw_electrodes(image);
    });


    for (let i = 0; i < document.getElementsByClassName("loadbutton").length; i++) {
        let load_button = document.getElementsByClassName("loadbutton")[i].id;
        document.getElementById(load_button).addEventListener("click", (_) => { loadJSONElectrodes(image, load_button); });
    }
}

function setupAdvancedButtons(image) {
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

    document.getElementById("JSON_SAVE").addEventListener("click", (_) => {
        const array = Array.from(image.electrode_pixels(), x => x);
        const file = new Blob([JSON.stringify({ map: array, state: state })], { type: 'text/plain' });
        let a = document.createElement("a");
        const url = URL.createObjectURL(file);
        a.href = url;
        a.download = "electrodes.json";
        document.body.appendChild(a);
        a.click();
        setTimeout(
            function () {
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            },
            0
        );
    });

    document.getElementById("JSON_LOAD").addEventListener("change", (event) => {
        const file = event.target.files[0];
        if (file) {
            const reader = new FileReader();
            reader.addEventListener("load", (e) => { loadJSONElectrodes(image, JSON.parse(reader.result)); });
            reader.readAsText(file);
        }
    });

    const view_button = document.getElementById('3d_view');
    view_button.addEventListener("click", (_) => {
        glob.refresh_canvas = 1;
        let electric_field_map = document.getElementById("electric_field_map");
        let new_efm = electric_field_map.cloneNode(false);
        electric_field_map.parentNode.replaceChild(new_efm, electric_field_map);
        if (!glob.three_dimension_view) {
            view_button.innerHTML = '2D</br>view';
            glob.gl = electric_field_map.getContext("webgl2");
            glob.three_dimension_view = 1;
        } else {
            view_button.innerHTML =  '3D</br>view';
            glob.gl = null;
            glob.three_dimension_view = 0;
        }
        setupElectricFieldMap(image);
        draw_electric_field(image);
    });

    document.getElementById("mz_ratio").addEventListener("input", (_) => { state.current_mz = parseFloat(document.getElementById("mz_ratio").value); });
    document.getElementById("ion_x_velocity").addEventListener("input", (_) => { state.current_x_velocity = parseFloat(document.getElementById("ion_x_velocity").value); });
    document.getElementById("ion_y_velocity").addEventListener("input", (_) => { state.current_y_velocity = parseFloat(document.getElementById("ion_y_velocity").value); });
    document.getElementById("randomizer").addEventListener("input", (_) => { state.random_ke = parseFloat(document.getElementById("randomizer").value); });
    document.getElementById("pulse_start").addEventListener("input", (_) => { state.pulse_starts[glob.active_electrode] = parseFloat(document.getElementById("pulse_start").value); });
    document.getElementById("pulse_end").addEventListener("input", (_) => { state.pulse_ends[glob.active_electrode] = parseFloat(document.getElementById("pulse_end").value); });
    document.getElementById("pressure").addEventListener("input", (_) => {
        state.pressure = parseFloat(document.getElementById("pressure").value);
        image.update_pressure(state.pressure);
    });
    document.getElementById("fdm_threshold").addEventListener("change", (_) => {
        state.fdm_accuracy = parseFloat(document.getElementById("fdm_threshold").value)
        image.update_fdm_threshold(state.fdm_accuracy);
    });  
    document.getElementById("bg_gas").addEventListener("change", (_) => {
        state.background_gas_mass = parseFloat(document.getElementById("bg_gas").value);
        image.update_gas_mass(state.background_gas_mass);
    });
    document.getElementById("sin_frequency").addEventListener("input", (a) => { state.electrode_frequencies[glob.active_electrode] = parseFloat(document.getElementById("sin_frequency").value); });
    document.getElementById("offset").addEventListener("input", (a) => {
        state.electrode_offsets[glob.active_electrode] = parseFloat(document.getElementById("offset").value);
        image.update_offsets(state.electrode_offsets);
    });
   
}


function setupPlayground() {
    let playground = document.getElementById("playground");
    let playgrounddiv = document.getElementById("playgrounddiv");
    let coll = document.getElementsByClassName("collapsible");
    for (let i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function () {
            this.classList.toggle("active");
            let content = this.nextElementSibling;
            if (content.style.display == "block") {
                content.style.display = "none";
            } else {
                content.style.display = "block";
            }
        });
        coll[i].nextElementSibling.style.display = "none";
    }
    playground.click();
    playground.click(); // nee
    playground.addEventListener("click", (event) => {
        if (playgrounddiv.style.backgroundColor == "rgb(243, 243, 243)") {
            if (glob.three_dimension_view) {
                document.getElementById("3d_view").click();
            }
            playgrounddiv.style.backgroundColor = "";
            playgrounddiv.style.boxShadow = "";
        } else {
            playgrounddiv.style.backgroundColor = "rgb(243, 243, 243)";
            playgrounddiv.style.boxShadow  = "0 2px 3px 0 #00000033, 0 2px 4px 0 #00000030";
        }
    })
    document.getElementById("advanced").click();
    document.getElementById("advanced").click();
    document.getElementById("sin_frequency").disabled = true;
    document.getElementById("offset").disabled = true;
    var header = document.getElementById("topnav");
    var sticky = playgrounddiv.offsetTop + header.offsetHeight;
    window.onscroll = () => {
        if (window.pageYOffset > sticky) {
            playgrounddiv.classList.add("sticky");
        } else {
            playgrounddiv.classList.remove("sticky");
        }
    };
    document.getElementById("dummyplayground").style.display = "none"
    document.getElementById("playgrounddiv").style.display = "inline-block"; // uncover the "playground" as soon as it's loaded
}


function refreshToggles() {
    document.getElementById("splattable").disabled = !Boolean(glob.active_electrode);
    document.getElementById("splattable").checked = state.splattable[glob.active_electrode];
    document.getElementById("dc_mode").disabled = !Boolean(glob.active_electrode);
    document.getElementById("dc_mode").checked = !state.dc_only[glob.active_electrode];
    document.getElementById("magnet_mode").checked = state.magnet_mode[glob.active_electrode];
    document.getElementById("sin_frequency").value = state.electrode_frequencies[glob.active_electrode];
    document.getElementById("offset").value = state.electrode_offsets[glob.active_electrode];
    document.getElementById("ion_x_velocity").value = state.current_x_velocity;
    document.getElementById("ion_y_velocity").value = state.current_y_velocity;
    document.getElementById("pressure").value = state.pressure;
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
    document.getElementById("voltage_input").value = ((glob.active_electrode == 0) ? 0 : state.electrode_volts[glob.active_electrode]);
}

function brushy_brush(image, map, width, height) {
    for (let x=0; x < width; x++)  {
        for (let y=0; y < height; y++) {
            image.brush(x, y, map[(x + width * y)]);
        }
    }
}

async function loadJSONElectrodes(image, imported) {
    let loaded_state;
    if (typeof imported === 'string' || imported instanceof String) {
        if (imported === 'clear') {
            return;
        } else if (!imported.endsWith("map")) {
            loaded_state = await import(`./maps/${imported}_map.json`).catch(console.error);
        } else {
            loaded_state = await import(`./maps/${imported}.json`).catch(console.error);
        }
    } else {
        loaded_state = imported;
    }
    let ls = loaded_state.state;
    brushy_brush(image, loaded_state.map, glob.map_width / 2, glob.map_height / 2); // load map
    Object.keys(state).forEach(k => { state[k] = ls[k]; }); // set global to new global
    glob.emap_altered = true;
    update_stuff(image);
    if (document.getElementById("playgrounddiv").style.backgroundColor == "") {
        document.getElementById("playground").click();
    }
    save_history(image);
}

function return_history(image) {
    if (!(glob.history1 == 0 && glob.history2 == 0)) {
        glob.history3 = glob.history2;
        glob.history2 = glob.history1;
        glob.history1 = 0;
        brushy_brush(image, glob.history3, glob.map_width / 2, glob.map_height / 2);
        glob.emap_altered = true;
    }
}

function update_stuff(image) {
    image.update_pressure(state.pressure);
    image.update_offsets(state.electrode_offsets);
    image.update_solids(state.splattable);
    image.update_magnets(state.magnet_mode);
    image.update_gas_mass(state.background_gas_mass);
    image.update_voltages(state.electrode_volts);
    refreshToggles();
}


export { setupButtons, setupAdvancedButtons, setupPlayground };