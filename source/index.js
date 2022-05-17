// We'll just make globals instead of passing a large state object. This keeps the code size a bit
// smaller and all important logic should be encapsulated in the wasm code.
import glob from './globals.js';
import { draw_electrodes, draw_over_electrodes, setupElectrodeMap } from './electrodes.js';
import { setupElectricFieldMap, draw_electric_field } from './electric_fields.js';
import { setupButtons } from './buttons.js';
import { width, height, max_time, scale, } from './constants.js';

const electrode_map = document.getElementById("electrode_map");
const electric_field_map = document.getElementById("electric_field_map");

export {
    electrode_map,
    electric_field_map,
};


function onTimerTick(image) {
    document.getElementById("electrode1").innerHTML = `${glob.electrode_volts[1]} V`;
    document.getElementById("electrode2").innerHTML = `${glob.electrode_volts[2]} V`;
    document.getElementById("electrode3").innerHTML = `${glob.electrode_volts[3]} V`;
    document.getElementById("electrode4").innerHTML = `${glob.electrode_volts[4]} V`;

    if (glob.mouse_on_electrode_map) {
        draw_over_electrodes(image);
    }
    if (glob.emap_altered) { // if the electrode pixels were changed, update the canvas
        glob.emap_altered = false;
        draw_electrodes(image);
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

function addStickyPlayground() {
    var header = document.getElementById("topnav");
    var playground = document.getElementById("playgrounddiv");
    var sticky = playground.offsetTop + header.offsetHeight;
    window.onscroll = () => {
        if (window.pageYOffset > sticky) {
            playground.classList.add("sticky");
        } else {
            playground.classList.remove("sticky");
        }
    };
}

function setupStuff() {
    var coll = document.getElementsByClassName("collapsible");
    for (let i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function () {
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.display == "block") {
                content.style.display = "none";
            } else {
                content.style.display = "block";
            }
        });
        coll[i].nextElementSibling.style.display = "none";
    }
    
    document.getElementById("playground").click();
    document.getElementById("playground").click(); // nee
    let playground = document.getElementById("playground");
    let playgrounddiv = document.getElementById("playgrounddiv");
    playground.addEventListener("click", (event) => {
        console.log(playgrounddiv.style.backgroundColor);
        if (playgrounddiv.style.backgroundColor == "rgb(243, 243, 243)") {
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
}

async function main() {
    document.getElementById("dummyplayground").style.display = "inline-block"
    const lib = await import("../pkg/index.js").catch(console.error);
    const image = new lib.Environment(
        width,
        height,
        scale,
        max_time,
        glob.electrode_volts,
        glob.splattable,
        glob.pressure,
        glob.background_gas_mass,
    );
    setupElectricFieldMap(image);
    setupElectrodeMap(image);
    setupButtons(image);
    setupStuff();
    draw_electrodes(image);
    draw_electric_field(image);
    highlightCurrentTopNavLink();
    addStickyPlayground();
    setInterval(onTimerTick.bind(null, image, electrode_map), 60);
    document.getElementById("dummyplayground").style.display = "none"
    document.getElementById("playgrounddiv").style.display = "inline-block"; // uncover the "playground" as soon as it's loaded
}

main();


