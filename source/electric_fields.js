import { electrode_map, electric_field_map } from './index.js';
import glob from './globals.js';

function drawFieldLines(image) {
    const context = electric_field_map.getContext("2d");
    const efpix = image.efpix();
    const width = image.width() * glob.cellSize;
    const height = image.height() * glob.cellSize;
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
            let sx = (x * glob.cellSize);
            let sy = (y * glob.cellSize);
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
    const pixels = image.efpix();
    const max = Math.max(...pixels);
    const min = Math.min(...pixels);
    for (let x = 0; x < glob.map_width; x++) {
        for (let y = 0; y < glob.map_height; y++) {
            let value = pixels[((y * glob.map_width) + x)];
            let newvalue = Math.floor((value - min) / (max - min) * 255.0);
            context.fillStyle = glob.cmap[newvalue];
            context.fillRect(x * glob.cellSize, y * glob.cellSize, glob.cellSize, glob.cellSize);
        }
    }
}

function setupElectricFieldMap(image) {
    electric_field_map.addEventListener('click', (event) => {
        const rect = electric_field_map.getBoundingClientRect();
        glob.current_x_position = (event.clientX - rect.left) / glob.cellSize;
        glob.current_y_position = (event.clientY - rect.top) / glob.cellSize;
        document.getElementById("x-position").value = glob.current_x_position;
        document.getElementById("y-position").value = glob.map_height - glob.current_y_position;
        fly_ion(image, glob.current_x_position, glob.current_y_position);
    });
    electric_field_map.addEventListener('mousemove', (event) => {
        const rect = electric_field_map.getBoundingClientRect();
        let x = event.clientX - rect.left;
        let y = event.clientY - rect.top;
        x = Math.floor(x / (glob.cellSize));
        y = Math.floor(y / (glob.cellSize));
        if (x >= 0 && y >= 0) {
            let ef = image.getef(x, y);
            x = `${x}`.padStart(3, " ");
            y = `${image.height() - y - 1}`.padStart(3, " ");
            ef = `${ef.toFixed(2)}`.padStart(5, " ");
            document.getElementById("status").innerHTML = `x: ${x} y: ${y}\nvolts: ${ef}`;
        }
    })
}


function fly_ion(image, x, y) {
    if (glob.electric_fields_calculated) {
        image.update_voltages(glob.electrode_volts);
        image.update_offsets(glob.electrode_offsets);
        glob.no_lines = false;
        let ix = x + 0.001;
        let iy = y + 0.001;
        let ivx = glob.current_x_velocity + 0.001;
        let ivy = -(glob.current_y_velocity + 0.001);
        let mz = glob.current_mz;
        glob.latest_ion = [ix, iy, ivx, ivy, mz, Math.random()];
        let ion_path = image.fly_ion(glob.latest_ion);
        glob.latest_ion_path = ion_path.slice(0, -1);
        let l = glob.latest_ion_path.length;
        plotLatestIonPath(electric_field_map);
        plotLatestIonPath(electrode_map);
        let fx = glob.latest_ion_path[l - 2].toFixed(0);
        let fy = image.height() - glob.latest_ion_path[l - 1].toFixed(0) - 1;
        let head = ((l / 3) >= glob.max_time) ? "Time!!": "Splat!";
        let position = `(${ix.toFixed(0)},${iy.toFixed(0)})->(${fx},${fy})`
        let collisions = ion_path[ion_path.length - 1];
        let splat = `${head} ${l / 3}ns ${mz}mz ${position} ${collisions} collisions`;
        glob.splatlog.push(splat);
        document.getElementById("splatlog").innerHTML = glob.splatlog.join('\n');
    } else {
        image.update_voltages(glob.electrode_volts);
        image.generate_electric_fields();
        glob.electric_fields_calculated = true;
        draw_electric_field(image);
    }
}


function plotLatestIonPath(canvas) {
    const context = canvas.getContext("2d");
    context.beginPath();
    context.strokeStyle = "#000";
    context.moveTo(glob.latest_ion_path[1] * glob.cellSize, glob.latest_ion_path[2] * glob.cellSize);
    var prev_x = 800;
    var prev_y = 800;
    for (let coordset = 0; coordset < glob.latest_ion_path.length; coordset += 3) {
        let x = glob.latest_ion_path[coordset + 1] * glob.cellSize;
        let y = glob.latest_ion_path[coordset + 2] * glob.cellSize;
        if ((Math.abs(x - prev_x) + Math.abs(y - prev_y)) > 0.25) {
            context.lineTo(x, y);
            prev_x = x;
            prev_y = y;
        }
    }
    context.stroke();
}

export { setupElectricFieldMap, drawFieldLines, draw_electric_field, fly_ion };