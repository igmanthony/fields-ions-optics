import { electrode_map } from './index.js';
import glob from './globals.js';
import { cellSize, max_time } from './constants.js';

function drawFieldLines(image) {
    const context = electric_field_map.getContext("2d");
    const efpix = image.efpix();
    const width = image.width() * cellSize;
    const height = image.height() * cellSize;
    for (let x = 1; x < width - 5; x += 5) {
        for (let y = 1; y < height - 5; y += 5) {
            let index = (y * (width / 3)) + x;
            context.beginPath();
            let t = efpix[index + 1];
            let b = efpix[index - 1];
            let r = efpix[index + width];
            let l = efpix[index - width];
            let hor = b - t;
            let ver = l - r;
            let hyp = Math.sqrt(hor * hor + ver * ver) / 8;
            let h = hor / hyp;
            let v = ver / hyp;
            context.strokeStyle = "#000";
            let sx = (x * cellSize);
            let sy = (y * cellSize);
            let ex = sx + h;
            let ey = sy + v;
            context.moveTo(sx, sy);
            context.lineTo(ex, ey);
            context.stroke();
            context.fill
        }
    }
}

function draw_electric_field(image) {
    const context = electric_field_map.getContext("2d");
    const pixels = image.electric_field_pixels();
    const width = image.width();
    const height = image.height();
    const size = cellSize;
    for (let x = 0; x < width; x++) {
        for (let y = 0; y < height; y++) {
            context.fillStyle = glob.cmap[pixels[((y * width) + x)]];
            context.fillRect(x * size, y * size, size, size);
        }
    }
}

function setupElectricFieldMap(image) {
    electric_field_map.addEventListener('click', (event) => {
        if (glob.electric_fields_calculated) {
            glob.no_lines = false;
            const rect = electric_field_map.getBoundingClientRect();
            let x = (event.clientX - rect.left) / cellSize;
            let y = (event.clientY - rect.top) / cellSize;
            let ix = x + 0.001;
            let iy = y + 0.001;
            let ivx = glob.current_x_velocity + 0.001;
            let ivy = -(glob.current_y_velocity + 0.001);
            let mz = glob.current_mz;
            glob.latest_ion = [ix, iy, ivx, ivy, mz, Math.random()];
            glob.latest_ion_path = image.fly_ion(glob.latest_ion);
            let l = glob.latest_ion_path.length;
            plotLatestIonPath(electric_field_map);
            plotLatestIonPath(electrode_map);
            ix = ix.toFixed(0);
            iy = iy.toFixed(0);
            ivy = ivy.toFixed(1);
            ivx = ivx.toFixed(1);
            let head = "Splat!"
            if ((l / 3) >= max_time) {
                head = "Time!!";
            }
            // let ft = (glob.latest_ion_path[l - 3] * 1000).toFixed(1);
            let fx = glob.latest_ion_path[l - 2].toFixed(0);
            let fy = image.height() - glob.latest_ion_path[l - 1].toFixed(0) - 1;
            let splat = `${head} ${l / 3 - 1} ns ${mz}mz (${ix},${iy})->(${fx},${fy})`;
            glob.splatlog.push(splat);
            document.getElementById("splatlog").innerHTML = glob.splatlog.join('\n');

        } else {
            image.update_voltages(glob.electrode_volts);
            console.log(`Updated voltages ${glob.electrode_volts}`);
            image.generate_electric_fields();
            console.log("Updated electrodes & generated electric fields");
            glob.electric_fields_calculated = true;
            draw_electric_field(image);
        }
    });
    electric_field_map.addEventListener('mousemove', (event) => {
        const rect = electric_field_map.getBoundingClientRect();
        let x = event.clientX - rect.left;
        let y = event.clientY - rect.top;
        x = Math.floor(x / (cellSize));
        y = Math.floor(y / (cellSize));
        if (x >= 0 && y >= 0) {
            let ef = image.getef(x, y);
            x = `${x}`.padStart(3, " ");
            y = `${image.height() - y - 1}`.padStart(3, " ");
            ef = `${ef.toFixed(2)}`.padStart(5, " ");
            document.getElementById("status").innerHTML = `x: ${x} y: ${y}\nvolts: ${ef}`;
        }
    })
}


function plotLatestIonPath(canvas) {
    const context = canvas.getContext("2d");
    context.beginPath();
    context.strokeStyle = "#000";
    context.moveTo(glob.latest_ion_path[1] * 3, glob.latest_ion_path[2] * 3);
    var prev_x = 800;
    var prev_y = 800;
    for (let coordset = 0; coordset < glob.latest_ion_path.length; coordset += 3) {
        let x = glob.latest_ion_path[coordset + 1] * 3;
        let y = glob.latest_ion_path[coordset + 2] * 3;
        if ((Math.abs(x - prev_x) + Math.abs(y - prev_y)) > 0.25) {
            context.lineTo(x, y);
            prev_x = x;
            prev_y = y;
        }
    }
    context.stroke();
}

export { setupElectricFieldMap, drawFieldLines, draw_electric_field };