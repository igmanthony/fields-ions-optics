import { electrode_map } from './index.js';
import glob from './globals.js';
import { electrode_colors, electrode_colors_100, cellSize } from './constants.js';

function draw_electrodes(image) {
    const scale = image.scale();
    const width = image.width() / scale;
    const height = image.height() / scale;
    const size = cellSize * scale;
    const pixels = image.electrode_pixels();
    const context = electrode_map.getContext("2d");
    // context.clearRect(0, 0, width * size, height * size);
    for (let x = 0; x < width; x++) {
        for (let y = 0; y < height; y++) {
            const pix = pixels[((y * width) + x)];
            if (glob.splattable[pix] == 1 && pix > 0) {
                context.fillStyle = electrode_colors_100[pix];
                context.fillRect(x * size, y * size, size, size);
                context.fillStyle = electrode_colors[pix];
                context.fillRect(x * size + 2, y * size + 2, size - 4, size - 4);
            } else {
                context.fillStyle = electrode_colors[pix];
                context.fillRect(x * size, y * size, size, size);
            }
        }
    }
}

function draw_over_electrodes(image) {
    if (!glob.dragging) {
        const scale = image.scale();
        const size = cellSize * scale;
        const rect = electrode_map.getBoundingClientRect();
        let x = glob.mouse_electrode_map_x - rect.left;
        let y = glob.mouse_electrode_map_y - rect.top;
        x = Math.floor(x / (cellSize * 2));
        y = Math.floor(y / (cellSize * 2));
        if (glob.previous_mouse_electrode_position_x2 != x || glob.previous_mouse_electrode_position_y2 != y) {
            // glob.emap_altered = true;
            draw_electrodes(image);
            glob.previous_mouse_electrode_position_x2 = x;
            glob.previous_mouse_electrode_position_y2 = y;
            const context = electrode_map.getContext("2d");
            context.fillStyle = electrode_colors[glob.active_electrode];
            if (glob.brush == 1) {
                context.fillRect(x * size, y * size, size * glob.brush, size * glob.brush);
            } else if (glob.brush == 2) {
                if (glob.initial_click[0] == -1 && glob.initial_click[1] == -1) {
                    context.fillRect(x * size, y * size, size, size);
                } else { // only works for mode 2 (straight lines)
                    let x0 = Math.floor(glob.initial_click[0] / (cellSize * 2));
                    let y0 = Math.floor(glob.initial_click[1] / (cellSize * 2));
                    let line = make_line(x0, y0, x, y );
                    for (let point of line) {
                        context.fillRect(point[0] * size, point[1] * size, size, size);
                    }
                }
            } else if (glob.brush == 3) {
                if (glob.initial_click[0] == -1 && glob.initial_click[1] == -1) {
                    context.fillRect(x * size, y * size, size, size);
                } else { // ellipse
                    let x0 = Math.floor(glob.initial_click[0] / (cellSize * 2));
                    let y0 = Math.floor(glob.initial_click[1] / (cellSize * 2));
                    let ellipse = make_ellipse(x0, y0, x, y );
                    for (let point of ellipse) {
                        context.fillRect(point[0] * size, point[1] * size, size, size);
                    }
                }
            } else if (glob.brush == 4) {
                if (glob.initial_click[0] == -1 && glob.initial_click[1] == -1) {
                    context.fillRect(x * size, y * size, size, size);
                } else { // only works for mode 2 (straight lines)
                    let x0 = Math.floor(glob.initial_click[0] / (cellSize * 2));
                    let y0 = Math.floor(glob.initial_click[1] / (cellSize * 2));
                    let box = make_box(x0, y0, x, y );
                    for (let point of box) {
                        context.fillRect(point[0] * size, point[1] * size, size, size);
                    }
                }
            }
            console.log("altered");
        }
    }
}

function make_line(p0x, p0y, p1x, p1y) {
    let points = [];
    let dx = p1x - p0x, dy = p1y - p0y;
    let N = Math.max(Math.abs(dx), Math.abs(dy));
    for (let step = 0; step <= N; step++) {
        let t = N === 0? 0.0 : step / N;
        let x = p0x + t * (p1x - p0x);
        let y = p0y + t * (p1y - p0y);
        points.push([Math.round(x), Math.round(y)]);
    }
    return points;
}

function make_box(p0x, p0y, p1x, p1y) {
    let top = make_line(p0x, p0y, p1x, p0y);
    let bottom = make_line(p0x, p1y, p1x, p1y);
    let left = make_line(p0x, p0y, p0x, p1y);
    let right = make_line(p1x, p0y, p1x, p1y);
    return [...top, ...left, ...right, ...bottom];
}

function make_ellipse(p0x, p0y, p1x, p1y) {
    let p2x = Math.abs(p0x - p1x) + p0x;
    let p2y = Math.abs(p0y - p1y) + p0y;
    let xs = Math.abs(p0x - p1x) * 2;
    let ys = Math.abs(p0y - p1y) * 2;
    let a = Math.abs(p2x - p0x);
    let b = Math.abs(p2y - p0y);
    let points = [];
    console.log(xs, ys, p0x, p0y, p2x, p2y, a, b);
    for (let y = p0y - p2y; y < p0y + p2y; y++) {
        let rows = [];
        for (let x = p0x - p2x; x < p0x + p2x; x++) {
            if (((x-p0x)*(x-p0x)/(a*a))+((y-p0y)*(y-p0y)/(b*b)) <= 1.0) {
                if (((x-p0x)*(x-p0x)/((a-1)*(a -1)))+((y-p0y)*(y-p0y)/((b-1)*(b-1))) > 1.0) {
                    points.push([Math.round(x), Math.round(y)]);
                }
            }
        }

    }
   return points;
}


function applyBrush(event, image) {
    if (!glob.dragging && event.type != 'click') { return; }

    glob.no_lines = true;
    const rect = electrode_map.getBoundingClientRect();
    let x = glob.mouse_electrode_map_x - rect.left;
    let y = glob.mouse_electrode_map_y - rect.top;
    glob.emap_altered = true;
    x = Math.floor(x / (cellSize * 2));
    y = Math.floor(y / (cellSize * 2));
    console.log("here", x, y);
    if ((x < 0) || (y < 0)) {
        return;
    }
    if (glob.brush == 1) {
        if (glob.previous_mouse_electrode_position_x != x || glob.previous_mouse_electrode_position_y != y) {
            console.log("applying")
            image.brush(x, y, glob.active_electrode);
            glob.previous_mouse_electrode_position_x = x;
            glob.previous_mouse_electrode_position_y = y;
        }
    } else if (glob.brush == 2 && !glob.dragging) {
        let x0 = Math.floor(glob.initial_click[0] / (cellSize * 2));
        let y0 = Math.floor(glob.initial_click[1] / (cellSize * 2));
        let line = make_line(x0, y0, x, y );
        for (let point of line) {
            image.brush(point[0], point[1], glob.active_electrode);
        }
    } else if (glob.brush == 3 && !glob.dragging) {
        let x0 = Math.floor(glob.initial_click[0] / (cellSize * 2));
        let y0 = Math.floor(glob.initial_click[1] / (cellSize * 2));
        let ellipse = make_ellipse(x0, y0, x, y);
        for (let point of ellipse) {
            image.brush(point[0], point[1], glob.active_electrode);
        }
    } else if (glob.brush == 4 && !glob.dragging) {
        let x0 = Math.floor(glob.initial_click[0] / (cellSize * 2));
        let y0 = Math.floor(glob.initial_click[1] / (cellSize * 2));
        let box = make_box(x0, y0, x, y );
        for (let point of box) {
            image.brush(point[0], point[1], glob.active_electrode);
        }
    }
}

function setupElectrodeMap(image) {
    electrode_map.addEventListener("mousedown", event => {
        glob.emap_altered = true;
        glob.dragging = true;
    });
    electrode_map.addEventListener("mouseup", event => {
        glob.emap_altered = true;
        glob.dragging = false;
    });
    electrode_map.addEventListener('click', (event) => {
        if (glob.brush == 1) {
            console.log("click apply", event.type);
            applyBrush(event, image);
        } else {
            // first click of line/square/circle
            if (glob.initial_click[0] == -1 && glob.initial_click[1] == -1) {
                const rect = electrode_map.getBoundingClientRect();
                glob.initial_click[0] = event.clientX - rect.left;
                glob.initial_click[1] = event.clientY - rect.top;

            } else { // second click of line/square/ circle
                applyBrush(event, image);
                glob.initial_click[0] = -1; // reset initial click status
                glob.initial_click[1] = -1; // reset initial click status
            }
        }
    });
    electrode_map.addEventListener('mousemove', (event) => {
        const rect = electrode_map.getBoundingClientRect();
        let x = Math.floor((event.clientX - rect.left) / (cellSize));
        let y = image.height() - (Math.floor((event.clientY - rect.top) / (cellSize))) - 1;
        document.getElementById("status").innerHTML = `x: ${x} y: ${y}`;
        glob.mouse_electrode_map_x = event.clientX;
        glob.mouse_electrode_map_y = event.clientY;
        applyBrush(event, image);
    });
    electrode_map.addEventListener('mouseover', (_) => { glob.mouse_on_electrode_map = true; });
    electrode_map.addEventListener('mouseout', (_) => {
        glob.mouse_on_electrode_map = false; 
        glob.dragging = false;
        glob.initial_click = [-1, -1],
        draw_electrodes(image);
    });
}

function makeTimeline(length) {
    let timeline = [];
    console.log("supposed timeline_length", length);
    for (let i=0; i < glob.electrode_volts.length; i++) {
        let dc_mode = glob.dc_only[i];
        console.log("dc_mode", dc_mode );
        let frequency = (glob.electrode_frequencies[i] * 1000) * 3.14159 * 2;
        let pulse_start = glob.pulse_starts[i];
        let pulse_end = glob.pulse_ends[i];
        for (let j = 0; j < length; j++) {
            if ((j <= pulse_end) && (j >= pulse_start)) {
                if (dc_mode) {
                    timeline.push(1.0);
                } else {
                    timeline.push( Math.sin( frequency * j * 1e-9) ); // 1 nm per step
                }
            } else {
                timeline.push(0.0);
            }
        }
    }
    let sum = 0;
    for (let i = 0; i < timeline.length; i++) {
        sum += timeline[i];
    }
    console.log("timeline made length: ", timeline.length, sum);
    return timeline;
}


export {
    draw_electrodes,
    draw_over_electrodes,
    setupElectrodeMap,
    makeTimeline,
};


