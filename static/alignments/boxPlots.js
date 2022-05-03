// Workaround for coloring individual boxpoints
// https://github.com/plotly/plotly.js/issues/2308

const zip = rows=>rows[0].map((_,c)=>rows.map(row=>row[c]))
function round(x, precision = 0) {
  let factor = Math.pow(10, precision);
  return Math.round(x * factor) / factor;
}

// Generate test data
let means = [5.25, 5.5, 6, 6.2, 6.6, 6.80] // average values
let sds = [3, 2, 2.5, 1, 8, 9] // standard deviations
// Map tuple(array) to measurement object
const Measurement = ([mean, sd]) => { return {mean, sd} } 
// Generate series from names, offset means by iteration index
const Series = (name, i) => {
    const offset = x => x + i
    return {
        name, 
        measurements: zip([means.map(offset), sds]).map(Measurement)
    }}
                                 
let series = ['A', 'B', 'C', 'D', 'E',].map(Series)

console.log(series)

const mean = (m) => m.mean
const sd = (m) => m.sd

// let ys = [].concat()
let seriesNames = series.map(s => s.name)

// Template boxplot with points we are trying to emulate
function traceBoxplotAll(ss) {
    return {
        name: 'all mean',
        x: [].concat(...ss.map(s => s.measurements.map(mean))),
        y: [].concat(...ss.map(
            s => Array(s.measurements.length).fill(s.name))),
        // y: ys(s), //series.indexOf(s),
        type: 'box', jitter: 1, pointpos: 2, boxpoints: 'all',
        orientation: 'h',
        marker: {color: 'rgba(7, 40, 89, .5)'},
        fillcolor: 'rgba(0, 0, 0, 0)', // transparent
        line: {color: 'rgba(0, 0, 0, 0.2)', },
    }
}

// Transparent fill over the interquartile range, to show we are painting
// over the template boxplot
function traceBoxplotIQR(ss) {
    return {
        name: 'means IQR',
        x: [].concat(...ss.map(s => s.measurements.map(mean))),
        y: [].concat(...ss.map(
            s => Array(s.measurements.length).fill(s.name))),
        // y: ys(s), //series.indexOf(s),
        type: 'box', boxpoints: false, orientation: 'h',
        fillcolor: 'rgba(7, 40, 89, .3)', 
        line: {width: 0}, // transparent
    }
}

const ys = (s) => {
    var y = Array(s.measurements.length)
    let i = series.indexOf(s)
    y.fill(i)
    return y
}

// repeatable pseudorandom generator from
// https://github.com/plotly/plotly.js/blob/master/src/traces/box/plot.js
var randSeed = 2000000000;
function rand() {
    var lastVal = randSeed;
    randSeed = (69069 * randSeed + 1) % 4294967296;
    // don't let consecutive vals be too close together
    // gets away from really trying to be random, in favor of better local uniformity
    if(Math.abs(randSeed - lastVal) < 429496729) return rand();
    return randSeed / 4294967296;
}

// jitter between 0.3 and 0.7 to stay in lane between boxplots w/o overlap
const jitter = (x) => x + 0.245 + round(rand() * 0.491, 3)
// Colorcoded point
function traceBoxpointsScatter(s) {
    return {
        name: s.name,
        x: s.measurements.map(mean),
        y: ys(s).map(jitter),
        yaxis: 'y2',
        type: 'scattergl', 
        mode: 'markers',
        opacity: 0.5,
        marker: {
            colorscale: 'Portland', cmin: 1, cmax: 10,
            color: s.measurements.map(sd),
        },
        text: s.measurements.map(m => `sd: ${m.sd}`),
        showlegend: false,
    }
}

let data = [
    traceBoxplotAll(series), 
    traceBoxplotIQR(series),
    ...series.map(traceBoxpointsScatter)
]
console.log(data)
let range = [10, -1.1]
let chart = document.getElementById('myDiv')

Plotly.newPlot( 
    chart,
    data, 
    {
          title: 'Box Plot with Colored Markers', 
        // showlegend: false,
        yaxis: {
            boxmode: 'grouped',
            range,
        },
        yaxis2: { 
            overlaying: 'y', 
            range,
            tickvals: [...Array(series.length).keys()],
            ticktext: series.map(s => s.name),
            side: 'right',
            zeroline: false,
            showgrid: false, // disable this to see axis going out-of-sync
            // without explicit re-synchronization in afterplot event handler
        }
    }
);

// Without this the Autoscale button gets the y and y2 out of sync.
chart.on('plotly_afterplot', synchronizeAxisRange)

const isEqual = (a1, a2) =>
  a1.length===a2.length && a1.every((v,i)=> v === a2[i])

function synchronizeAxisRange() {
    let l = chart.layout
    let r = l.yaxis.range.map(r => round(r, 3))
    const setRange = (axis, r) => { axis.autorange = false; axis.range = r }
    if (!isEqual(r, l.yaxis2.range)) {
        setRange(l.yaxis, r)
        setRange(l.yaxis2, r)
        Plotly.relayout(chart, l)
  }
}