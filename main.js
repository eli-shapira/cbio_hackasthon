var request = require('request');
var { JSDOM } = require('jsdom');
var jQuery = require('jquery');


function parse_uniport(res, prot){
    let sum = 'protein:' + prot + '\n';
    sum += res + '\n' + '##################################################' + '\n';
    return sum;
}


const filters = {
    returnALL: _ => true
}


// var human_list = require('./human_prot_list');
// var yeast_list = require('./yeast_prot_list');
var dolphin_list = require('./dolphin_prot_list');
var jelly_list = require('./jelly_prot_list');


const UNIPORT_PATH = prot => `https://www.uniprot.org/uniprot/${prot}.txt`;


function getHTML(base_url, acc) {
    let url = base_url(acc);
    return new Promise((resolve, reject) => {
        request({url: url}, function (error, response, body) {
            if (error) reject(error);
            else if (parseInt(response.statusCode / 100) !== 2) {
                reject(response.statusMessage);
            }
            else resolve([body, acc]);
        });
    })
}


function initJQuery(responseBody) {
    let dom = new JSDOM(responseBody);
    return jQuery(dom.window);
}


function sub_run(base_url, list, parse_func){
    let promisesArr = list.map( route => {
        return new Promise((res, rej) => {
            getHTML(base_url, route).then(response => {
                res(parse_func(response[0], response[1]));
            }).catch(rej);
        });
    });
    // console.log(promisesArr);
    return Promise.all(promisesArr);
}


async function run_all(base_url, list, sliceParam, parse_func, filter_results) {
    /**
     * @param {Array} list - list of params to get access in GEO
     * @param {Int} sliceParam - how many samples to run at each batch
     * @param {Function} parse_func - one of parser module functions
     * @param {Function} filter_results - in the form (result_of_parse_func) => some_condition on it
     */
    let numOfSlices = parseInt(list.length / sliceParam);
    // listOfLists = []
    for (let i = 0; i <= numOfSlices; i++) {
    // for (let i = 0; i <= 2; i++) {
        // listOfLists.push();
        let sumList = await sub_run(base_url, list.slice(sliceParam * i, sliceParam * (i + 1)), parse_func);
        console.log(sumList.filter(filter_results).join('\n'));
    }
}


program_arguments = process.argv;


if (program_arguments.length > 2 && program_arguments[2] === 'stdin') {
    var stdin = process.openStdin();
    stdin.addListener("data", function(d) {
        var input = d.toString().trim();
        if (input.indexOf("^ ") === 0) { 
            getHTML(input.split(' ')[1]).then(response => {
                console.log(response.length);
                initJQuery(response);
            });
        }
        else if (input === 'exit') process.exit()
        console.log(input);
    });
} else {
    // run_all(UNIPORT_PATH, yeast_list, 10, parse_uniport, filters.returnALL).catch(err => console.log(err));
    // run_all(UNIPORT_PATH, human_list, 10, parse_uniport, filters.returnALL).catch(err => console.log(err));
    // run_all(UNIPORT_PATH, jelly_list, 10, parse_uniport, filters.returnALL).catch(err => console.log(err));
    run_all(UNIPORT_PATH, dolphin_list, 10, parse_uniport, filters.returnALL).catch(err => console.log(err));
    // run_all(UNIPORT_PATH, ['P02730'], 10, parse_uniport, filters.returnALL);
}
