import { glob, state } from './globals.js';
import { draw_electrodes, draw_over_electrodes, setupElectrodeMap } from './electrodes.js';
import { setupElectricFieldMap, draw_electric_field } from './fields.js';
import { setupAdvancedButtons, setupButtons, setupPlayground } from './buttons.js';


function onTimerTick(image) {
    if (!glob.dragging && glob.brush != 1) {
        glob.to_apply = []; // clear the points to apply to prevent oddities
    }
    for (let i = 1; i < 5; i++) {
        document.getElementById(`electrode${i}`).innerHTML = `${state.electrode_volts[i]} ${state.magnet_mode[i] ? "T" : "V"}`;
    }
    document.getElementById("set_voltage").innerHTML = state.magnet_mode[glob.active_electrode] ? "Set Teslas" : "Set Volts";
    if (glob.mouse_on_electrode_map) {
        draw_over_electrodes(image);
    }
    if (glob.emap_altered) { // if the electrode pixels were changed, update the canvas
        draw_electrodes(image);
        glob.emap_altered = false;
    }
    if (glob.voltages_were_updated) {
        draw_electric_field(image);
        glob.voltages_were_updated = false;
    }
}

function highlightCurrentTopNavLink() {
    let loc = location.pathname.split('/');
    let current = loc[loc.length - 1];
    if (current !== "") {
        let menuItems = document.querySelectorAll('.topnav a');
        for (let i = 0, len = menuItems.length; i < len; i++) {
            if (menuItems[i].getAttribute("href").indexOf(current) == 0) { // removed !== -1
                menuItems[i].className += "current_page";
            }
        }
    }
}

async function main() {
    document.getElementById("dummyplayground").style.display = "inline-block"
    const lib = await import("../pkg/index.js").catch(console.error);
    const image = new lib.Environment(
        glob.map_width,
        glob.map_height,
        glob.scale,
        state.max_time,
        state.electrode_volts,
        state.splattable,
        state.pressure,
        state.background_gas_mass,
    );
    let electrode_map = document.getElementById("electrode_map");
    let electric_field_map = document.getElementById("electric_field_map");
    if ((window.screen.height * window.devicePixelRatio < 1070)) {
        glob.cellSize = 2; // usually 3
        glob.canvas_height = 400; // ususally 600
        glob.canvas_width = 400; // ususally 600
        electric_field_map.width = '400';
        electric_field_map.height = '400';
        electrode_map.width = '400';
        electrode_map.height = '400';
        document.getElementById("status").style.bottom = "-445px";
        document.getElementById("status").style.left = "-380px";
        document.getElementById("splatlog").style.width = "280px";
    };
    highlightCurrentTopNavLink();
    setupPlayground();
    setupButtons(image);
    setupAdvancedButtons(image);
    setupElectricFieldMap(image);
    setupElectrodeMap(image);
    draw_electrodes(image);
    draw_electric_field(image);
    setInterval(onTimerTick.bind(null, image, electrode_map), 60);
}

main();


