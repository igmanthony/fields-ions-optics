import glob from './globals.js';
import {draw3d} from './render.js';


function draw_electric_field(image) {
    const pixels = image.efpix();
    glob.efmap_pixels = pixels;
    let electric_field_map = document.getElementById("electric_field_map");
    if (glob.three_dimension_view) {
        if (glob.gl === null) {
            glob.gl = electric_field_map.getContext("webgl2")
        }
        draw3d(image);
    } else {
        const context = electric_field_map.getContext("2d");
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
}

function setupElectricFieldMap(image) {
    let electric_field_map = document.getElementById("electric_field_map");
    glob.three_d_zoom_level = 0;
    if (glob.three_dimension_view) {
        electric_field_map.addEventListener('click', (event) => {
            const rect = electric_field_map.getBoundingClientRect();
            glob.current_x_position = (event.clientX - rect.left) / glob.cellSize;
            glob.current_y_position = (event.clientY - rect.top) / glob.cellSize;
            document.getElementById("x-position").value = glob.current_x_position;
            document.getElementById("y-position").value = glob.map_height - glob.current_y_position;
        });
        electric_field_map.addEventListener('mousemove', (event) => {
            if (glob.dragging) {
                let initial_y = glob.mouse_ef_map_y;
                let initial_x = glob.mouse_ef_map_x;
                glob.mouse_ef_map_y = event.clientX; // opposites?
                glob.mouse_ef_map_x = event.clientY; // opposites?
                glob.delta_ef_map_x += initial_x - glob.mouse_ef_map_x;
                glob.delta_ef_map_y += initial_y - glob.mouse_ef_map_y;
            }
        })

        electric_field_map.addEventListener("mousedown", event => {
            glob.dragging = true;
            glob.mouse_ef_map_y = event.clientX; // opposites?
            glob.mouse_ef_map_x = event.clientY; // opposites?
        });
        electric_field_map.addEventListener("mouseup", event => { glob.dragging = false; });
        electric_field_map.addEventListener('wheel', event => {
            event.preventDefault();
            glob.three_d_zoom_level += event.deltaY / 2.0;
        });
    } else {
        let electric_field_map = document.getElementById("electric_field_map");
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
        let seed = Math.random();
        if (glob.random_ke != 0) {
            let angle = seed * 2 * Math.PI;
            let ke = glob.random_ke * 1.60218E-19; // convert from μeV to joules
            let v = Math.sqrt((2 * ke) / (mz * 1.660538921E-27)) / 1000.0; // in m/ms
            ivx = ivx + (v * Math.cos(angle));
            ivy = ivy + (v * Math.sin(angle));
        }
        glob.latest_ion = [ix, iy, ivx, ivy, mz, seed];
        let ion_path = image.fly_ion(glob.latest_ion);
        glob.latest_ion_path = ion_path.slice(0, -1);
        let l = glob.latest_ion_path.length;
        if (glob.three_dimension_view) {
            plotLatestIonPath(glob.fake_canvas, 1);
            const gl = glob.gl;
            let texture = glob.fake_canvas.getImageData(0, 0, glob.map_width, glob.map_height);
            gl.bindTexture(gl.TEXTURE_2D, gl.createTexture());
            gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texture);
            gl.generateMipmap(gl.TEXTURE_2D);
        } else {
            const electric_field_map = document.getElementById("electric_field_map");
            const ctx = electric_field_map.getContext("2d");
            plotLatestIonPath(ctx, glob.cellSize);
        }
        let electrode_map = document.getElementById("electrode_map");
        const emc = electrode_map.getContext("2d");
        plotLatestIonPath(emc, glob.cellSize);
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


function plotLatestIonPath(context, scaleFactor) {
    context.beginPath();
    context.strokeStyle = "#000";
    context.moveTo(glob.latest_ion_path[1] * scaleFactor, glob.latest_ion_path[2] * scaleFactor);
    var prev_x = 800;
    var prev_y = 800;
    for (let coordset = 0; coordset < glob.latest_ion_path.length; coordset += 3) {
        let x = glob.latest_ion_path[coordset + 1] * scaleFactor;
        let y = glob.latest_ion_path[coordset + 2] * scaleFactor;
        if ((Math.abs(x - prev_x) + Math.abs(y - prev_y)) > 0.25) {
            context.lineTo(x, y);
            prev_x = x;
            prev_y = y;
        }
    }
    context.stroke();
}

export { setupElectricFieldMap, draw_electric_field, fly_ion };