const turbo = [
    "#30123b", "#311542", "#32184a", "#341b51", "#351e58", "#36215f", "#372365", "#38266c",
    "#392972", "#3a2c79", "#3b2f7f", "#3c3285", "#3c358b", "#3d3791", "#3e3a96", "#3f3d9c",
    "#4040a1", "#4043a6", "#4145ab", "#4148b0", "#424bb5", "#434eba", "#4350be", "#4353c2",
    "#4456c7", "#4458cb", "#455bce", "#455ed2", "#4560d6", "#4563d9", "#4666dd", "#4668e0",
    "#466be3", "#466de6", "#4670e8", "#4673eb", "#4675ed", "#4678f0", "#467af2", "#467df4",
    "#467ff6", "#4682f8", "#4584f9", "#4587fb", "#4589fc", "#448cfd", "#438efd", "#4291fe",
    "#4193fe", "#4096fe", "#3f98fe", "#3e9bfe", "#3c9dfd", "#3ba0fc", "#39a2fc", "#38a5fb",
    "#36a8f9", "#34aaf8", "#33acf6", "#31aff5", "#2fb1f3", "#2db4f1", "#2bb6ef", "#2ab9ed",
    "#28bbeb", "#26bde9", "#25c0e6", "#23c2e4", "#21c4e1", "#20c6df", "#1ec9dc", "#1dcbda",
    "#1ccdd7", "#1bcfd4", "#1ad1d2", "#19d3cf", "#18d5cc", "#18d7ca", "#17d9c7", "#17dac4",
    "#17dcc2", "#17debf", "#18e0bd", "#18e1ba", "#19e3b8", "#1ae4b6", "#1be5b4", "#1de7b1",
    "#1ee8af", "#20e9ac", "#22eba9", "#24eca6", "#27eda3", "#29eea0", "#2cef9d", "#2ff09a",
    "#32f197", "#35f394", "#38f491", "#3bf48d", "#3ff58a", "#42f687", "#46f783", "#4af880",
    "#4df97c", "#51f979", "#55fa76", "#59fb72", "#5dfb6f", "#61fc6c", "#65fc68", "#69fd65",
    "#6dfd62", "#71fd5f", "#74fe5c", "#78fe59", "#7cfe56", "#80fe53", "#84fe50", "#87fe4d",
    "#8bfe4b", "#8efe48", "#92fe46", "#95fe44", "#98fe42", "#9bfd40", "#9efd3e", "#a1fc3d",
    "#a4fc3b", "#a6fb3a", "#a9fb39", "#acfa37", "#aef937", "#b1f836", "#b3f835", "#b6f735",
    "#b9f534", "#bbf434", "#bef334", "#c0f233", "#c3f133", "#c5ef33", "#c8ee33", "#caed33",
    "#cdeb34", "#cfea34", "#d1e834", "#d4e735", "#d6e535", "#d8e335", "#dae236", "#dde036",
    "#dfde36", "#e1dc37", "#e3da37", "#e5d838", "#e7d738", "#e8d538", "#ead339", "#ecd139",
    "#edcf39", "#efcd39", "#f0cb3a", "#f2c83a", "#f3c63a", "#f4c43a", "#f6c23a", "#f7c039",
    "#f8be39", "#f9bc39", "#f9ba38", "#fab737", "#fbb537", "#fbb336", "#fcb035", "#fcae34",
    "#fdab33", "#fda932", "#fda631", "#fda330", "#fea12f", "#fe9e2e", "#fe9b2d", "#fe982c",
    "#fd952b", "#fd9229", "#fd8f28", "#fd8c27", "#fc8926", "#fc8624", "#fb8323", "#fb8022",
    "#fa7d20", "#fa7a1f", "#f9771e", "#f8741c", "#f7711b", "#f76e1a", "#f66b18", "#f56817",
    "#f46516", "#f36315", "#f26014", "#f15d13", "#ef5a11", "#ee5810", "#ed550f", "#ec520e",
    "#ea500d", "#e94d0d", "#e84b0c", "#e6490b", "#e5460a", "#e3440a", "#e24209", "#e04008",
    "#de3e08", "#dd3c07", "#db3a07", "#d93806", "#d73606", "#d63405", "#d43205", "#d23005",
    "#d02f04", "#ce2d04", "#cb2b03", "#c92903", "#c72803", "#c52602", "#c32402", "#c02302",
    "#be2102", "#bb1f01", "#b91e01", "#b61c01", "#b41b01", "#b11901", "#ae1801", "#ac1601",
    "#a91501", "#a61401", "#a31201", "#a01101", "#9d1001", "#9a0e01", "#970d01", "#940c01",
    "#910b01", "#8e0a01", "#8b0901", "#870801", "#840701", "#810602", "#7d0502", "#7a0402",
];

const electrode_colors = [ "#ffffff", "#f44336", "#4CAF50", "#2196F3", "#795548", ];

const electrode_colors_100 = [ "#ffffff", "#FFCDD2", "#C8E6C9", "#BBDEFB", "#D7CCC8", ];

const electrode_colors_300 = [ "#ffffff", "#E57373", "#81C784", "#64B5F6", "#A1887F", ];

const default_splats = [-1, 0, 0, 0, 0];
const default_pressure = 0;
const default_dcs = [ true, true, true, true, true, ];
const default_frequencies = [0, 100, 100, 100, 100];
const default_offs = [0, 0, 0, 0, 0];

const einzel_lens_values = [0, 70, 100, 1, -1];
const tof_values = [0, 1000, -500, 0, 1020];
const tof_splats = [-1, 0, 0, 1, 0];
const quadrupole_values = [0, 40, 1, -40, 1];
const quadrupole_frequencies = [0, 110, 110, 110, 110];
const quadrupole_offsets = [0, 10, 0, 10, 0];
const cyclotron_frequencies = [0, 153.55994, 0, 153.55994, 0];
const cyclotron_magnets = [0, 0, 1, 0, 0];
const funnel_values = [0, 150, 350, -150, -400];
const funnel_frequencies = [0, 70, 70, 70, 70];
const funnel_pressure = 30;
const quadrupole_dcs = [ true, false, true, false, true ];

const green200 = '#A5D6A7';
const amber200 = '#FFE082';
const blue700 = '#1976D2';

export { 
    turbo,
    electrode_colors,
    electrode_colors_300,
    electrode_colors_100,
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
    cyclotron_magnets,
    tof_values,
    tof_splats,
    funnel_values,
    funnel_frequencies,
    funnel_pressure,
    green200,
    amber200,
    blue700,
    cyclotron_frequencies,
};