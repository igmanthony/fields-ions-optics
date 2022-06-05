import glob from './globals.js';
import { electrode_colors, electrode_colors_100 } from './constants.js';

function draw_electrodes(image) {
    const scale = image.scale();
    const width = image.width() / scale;
    const height = image.height() / scale;
    const size = glob.cellSize * scale;
    const pixels = image.electrode_pixels();
    let electrode_map = document.getElementById("electrode_map");
    const context = electrode_map.getContext("2d");
    for (let x = 0; x < width; x++) {
        for (let y = 0; y < height; y++) {
            const pix = pixels[((y * width) + x)];
            if (glob.magnet_mode[pix] == 1 && pix > 0) {
                context.fillStyle = electrode_colors_100[pix];
                context.fillRect(x * size, y * size, size, size);
                if (x % 3 == 0 && y % 3 == 0) {
                    context.fillStyle = electrode_colors[pix];
                    if (glob.electrode_volts[pix] > 0) {
                        for (let i = 0; i < size; i++) {
                            context.fillRect(x * size + i, y * size + i, 1, 1);
                            context.fillRect(x * size - i + (size-1), y * size + i, 1, 1);
                        }
                    } else {
                        context.fillRect(x * size + 2, y * size + 2, 2, 2);
                    }
                }
            } else if (glob.splattable[pix] == 1 && pix > 0) {
                context.fillStyle = electrode_colors[0];
                context.fillRect(x * size, y * size, size, size);
                context.fillStyle = electrode_colors_100[pix];
                context.fillRect(x * size + 2, y * size + 2, 2, 2);
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
        const size = glob.cellSize * scale;
        let electrode_map = document.getElementById("electrode_map");
        const rect = electrode_map.getBoundingClientRect();
        let x = glob.mouse_electrode_map_x - rect.left;
        let y = glob.mouse_electrode_map_y - rect.top;
        x = Math.floor(x / (glob.cellSize * 2));
        y = Math.floor(y / (glob.cellSize * 2));
        if (glob.previous_mouse_electrode_position_x2 != x || glob.previous_mouse_electrode_position_y2 != y) {
            // glob.emap_altered = true;
            draw_electrodes(image);
            glob.previous_mouse_electrode_position_x2 = x;
            glob.previous_mouse_electrode_position_y2 = y;
            const context = electrode_map.getContext("2d");
            context.fillStyle = electrode_colors[glob.active_electrode];
            if (glob.brush == 1 || glob.brush == 5) {
                context.fillRect(x * size, y * size, size * 1, size * 1);
            } else if (glob.brush == 2) {
                if (glob.initial_click[0] == -1 && glob.initial_click[1] == -1) {
                    context.fillRect(x * size, y * size, size, size);
                } else { // only works for mode 2 (straight lines)
                    let x0 = Math.floor(glob.initial_click[0] / (glob.cellSize * 2));
                    let y0 = Math.floor(glob.initial_click[1] / (glob.cellSize * 2));
                    let line = make_line(x0, y0, x, y );
                    for (let point of line) {
                        context.fillRect(point[0] * size, point[1] * size, size, size);
                    }
                }
            } else if (glob.brush == 3) {
                if (glob.initial_click[0] == -1 && glob.initial_click[1] == -1) {
                    context.fillRect(x * size, y * size, size, size);
                } else { // ellipse
                    let x0 = Math.floor(glob.initial_click[0] / (glob.cellSize * 2));
                    let y0 = Math.floor(glob.initial_click[1] / (glob.cellSize * 2));
                    let ellipse = make_ellipse(x0, y0, x, y );
                    for (let point of ellipse) {
                        context.fillRect(point[0] * size, point[1] * size, size, size);
                    }
                }
            } else if (glob.brush == 4) {
                if (glob.initial_click[0] == -1 && glob.initial_click[1] == -1) {
                    context.fillRect(x * size, y * size, size, size);
                } else { // only works for mode 2 (straight lines)
                    let x0 = Math.floor(glob.initial_click[0] / (glob.cellSize * 2));
                    let y0 = Math.floor(glob.initial_click[1] / (glob.cellSize * 2));
                    let box = make_box(x0, y0, x, y );
                    for (let point of box) {
                        context.fillRect(point[0] * size, point[1] * size, size, size);
                    }
                }
            }
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

function floodfill(x0, y0, points, width, height) {
    let updated_points = [];
    let start_color = points[(y0 * width) + x0];
    let Q = [];
    Q.push([x0, y0]);
    while (Q.length != 0) {
        let n = Q.pop();
        let y = n[1];
        let x = n[0];
        let current_color = points[(y * width) + x];
        if (current_color == start_color) {
            updated_points.push((y * width) + x);
            for (let p of [[x + 1, y], [x - 1, y], [x, y + 1], [x, y - 1]]) {
                if (0 <= p[0] && p[0] <= width-1 && 0 <= p[1] && p[1] <= height - 1) {
                    if (!updated_points.includes((p[1] * width) + p[0] )) {
                        Q.push(p);
                    }
                }
            }
        }
    }
    return updated_points;
}

function save_history(image) {
    glob.history1 = glob.history2;
    glob.history2 = glob.history3;
    glob.history3 = image.electrode_pixels();
} 

function return_history(image) {
    if (glob.history1 == 0 && glob.history2 == 0) {
        return
    }
    glob.history3 = glob.history2;
    glob.history2 = glob.history1;
    glob.history1 = 0;
    let width = glob.map_width / 2;
    let height = glob.map_height / 2;
    for (let x=0; x < width; x++)  {
        for (let y=0; y < height; y++) {
            let color = glob.history3[(x + width * y)];
            image.brush(x, y, color);
        }
    }
    glob.emap_altered = true;
}

function applyBrush(event, image) {
    if (glob.dragging) {
        glob.to_apply.push([glob.mouse_electrode_map_x, glob.mouse_electrode_map_y]);
    } else if (event.type == 'click') {
        glob.no_lines = true;
        let electrode_map = document.getElementById("electrode_map");
        const rect = electrode_map.getBoundingClientRect();
        let x = glob.mouse_electrode_map_x - rect.left;
        let y = glob.mouse_electrode_map_y - rect.top;
        glob.emap_altered = true;
        x = Math.floor(x / (glob.cellSize * 2));
        y = Math.floor(y / (glob.cellSize * 2));
        if ((x < 0) || (y < 0)) {
            return;
        }
        if (glob.brush == 1) {
            if (glob.previous_mouse_electrode_position_x != x || glob.previous_mouse_electrode_position_y != y) {
                image.brush(x, y, glob.active_electrode);
                glob.previous_mouse_electrode_position_x = x;
                glob.previous_mouse_electrode_position_y = y;
            }
        } else if (glob.brush == 2 && !glob.dragging) {
            let x0 = Math.floor(glob.initial_click[0] / (glob.cellSize * 2));
            let y0 = Math.floor(glob.initial_click[1] / (glob.cellSize * 2));
            let line = make_line(x0, y0, x, y );
            for (let point of line) {
                image.brush(point[0], point[1], glob.active_electrode);
            }
        } else if (glob.brush == 3 && !glob.dragging) {
            let x0 = Math.floor(glob.initial_click[0] / (glob.cellSize * 2));
            let y0 = Math.floor(glob.initial_click[1] / (glob.cellSize * 2));
            let ellipse = make_ellipse(x0, y0, x, y);
            for (let point of ellipse) {
                image.brush(point[0], point[1], glob.active_electrode);
            }
        } else if (glob.brush == 4 && !glob.dragging) {
            let x0 = Math.floor(glob.initial_click[0] / (glob.cellSize * 2));
            let y0 = Math.floor(glob.initial_click[1] / (glob.cellSize * 2));
            let box = make_box(x0, y0, x, y );
            for (let point of box) {
                image.brush(point[0], point[1], glob.active_electrode);
            }
        } else if (glob.brush == 5 && !glob.dragging) {
            const scale = image.scale();
            const width = image.width() / scale;
            const height = image.height() / scale;
            let points = image.electrode_pixels();
            let fill = floodfill(x, y, points, width, height);
            for (let point of fill) {
                image.brush(point % width, Math.floor(point / width), glob.active_electrode);
            }
        }
        save_history(image);
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
        if (glob.brush == 1) {
            const rect = electrode_map.getBoundingClientRect();
            for (let p of glob.to_apply) {
                let x = p[0] - rect.left;
                let y = p[1] - rect.top;
                x = Math.round(x / (glob.cellSize * 2));
                y = Math.round(y / (glob.cellSize * 2));
                image.brush(x, y, glob.active_electrode);
            }
            glob.to_apply = [];
        }
    });
    electrode_map.addEventListener('click', (event) => {
        if (glob.brush == 1 || glob.brush == 5) {
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
        const size = glob.cellSize * glob.scale;
        glob.rect = [rect.left, rect.top];
        let x = Math.floor((event.clientX - rect.left) / (glob.cellSize));
        let y = image.height() - (Math.floor((event.clientY - rect.top) / (glob.cellSize))) - 1;
        document.getElementById("status").innerHTML = `x: ${x} y: ${y}`;
        glob.mouse_electrode_map_x = event.clientX;
        glob.mouse_electrode_map_y = event.clientY;
        applyBrush(event, image);
        if (glob.dragging && glob.brush == 1) {
            const context = electrode_map.getContext("2d");
            context.fillStyle = electrode_colors[glob.active_electrode];
            for (let p of glob.to_apply) {
                let x = p[0] - rect.left;
                let y = p[1] - rect.top;
                x = Math.round(x / (glob.cellSize * 2));
                y = Math.round(y / (glob.cellSize * 2));
                context.fillRect(x * size, y * size, size, size);
            }
        }
    });
    electrode_map.addEventListener('mouseover', (_) => { glob.mouse_on_electrode_map = true; });
    electrode_map.addEventListener('mouseout', (_) => {
        glob.mouse_on_electrode_map = false; 
        glob.dragging = false;
        glob.initial_click = [-1, -1];
        if (glob.brush == 1) {
            const rect = electrode_map.getBoundingClientRect();
            for (let p of glob.to_apply) {
                let x = p[0] - rect.left;
                let y = p[1] - rect.top;
                x = Math.round(x / (glob.cellSize * 2));
                y = Math.round(y / (glob.cellSize * 2));
                image.brush(x, y, glob.active_electrode);
            }
            glob.to_apply = [];
        }
        draw_electrodes(image);
    });
}

export {
    draw_electrodes,
    draw_over_electrodes,
    setupElectrodeMap,
    return_history,
    save_history,
};


