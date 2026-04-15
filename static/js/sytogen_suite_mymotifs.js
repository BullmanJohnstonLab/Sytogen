// -----------------------------------------------------------------------------
// MyMotifs
// -----------------------------------------------------------------------------

function get_my_motif_choice(this_obj) {

    var choice = $('input[name="optradio"]:checked').val();

    // Modify your REBASE file
    // Create your REBASE file

    $('#motifs_dict').html(JSON.stringify({}));

    if (choice == 'Create your REBASE file') {
        $('#upload_rebase_button').remove();
        $('#upload_input_file').remove();
        $('#input_body').append('<div id="rebase_placeholder"></div>');
        $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
        $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="save_rebase()">Save rebase</button>');
    } else {
        reload()
    }

}

function draw_table_of_motifs(n_dict) {

    if ($("#rebase_output_table").length == 0) {

        $('#input_body').append(
            '<table class="table" id="rebase_output_table">' +
            '<thead>' +
            '<tr>' +
            '<th scope="col">Type</th>' +
            '<th scope="col">Recognition motif</th>' +
            '<th scope="col">Methylated base (+)</th>' +
            '<th scope="col">Methylated base (+) type</th>' +
            '<th scope="col">Methylated base (-)</th>' +
            '<th scope="col">Methylated base (-) type</th>' +
            '<th scope="col">Available</th>' +
            '<th scope="col">Confirm motif</th>' +
            '</tr>' +
            '</thead>' +
            '<tbody>' +
            '</tbody>' +
            '</table>'
        )

    }



    var table = document.getElementById('rebase_output_table');

    for (key1 of Object.keys(n_dict)) {
        
        if ((Object.keys(n_dict[key1]).includes('meth_base')) & (!Object.keys(n_dict[key1]).includes('comp_meth_base'))) {
            n_dict[key1]['comp_meth_base'] = "-"
            n_dict[key1]['comp_meth_type'] = "-"
        } else if ((!Object.keys(n_dict[key1]).includes('meth_base')) & (Object.keys(n_dict[key1]).includes('comp_meth_base'))) {
            n_dict[key1]['meth_base'] = "-"
            n_dict[key1]['meth_type'] = "-"
        }

        console.log("Here: ", n_dict[key1]);

        if (!Object.keys(n_dict[key1]).includes('enz_type')) {
            n_dict[key1]['enz_type'] = "-1"
        }

        var rec_seq_repr = {};

        for (let i = 0; i < n_dict[key1]['rec_seq'].length; i++) {

            rec_seq_repr[i] = n_dict[key1]['rec_seq'].charAt(i);
        }

        if (!Object.keys(n_dict[key1]).includes('meth_base')) {
            n_dict[key1]['meth_base'] = "-"
        } else {
            if (n_dict[key1]['meth_base'] != "-") {
                rec_seq_repr[Number(n_dict[key1]['meth_base']) - 1] = rec_seq_repr[Number(n_dict[key1]['meth_base']) - 1].bold();
            }

        }

        if (!Object.keys(n_dict[key1]).includes('meth_type')) {
            n_dict[key1]['meth_type'] = "-"
        }

        if (!Object.keys(n_dict[key1]).includes('comp_meth_base')) {
            n_dict[key1]['comp_meth_base'] = "-"
        } else {
            if (n_dict[key1]['comp_meth_base'] != "-") {
                rec_seq_repr[Number(n_dict[key1]['comp_meth_base']) - 1] = rec_seq_repr[Number(n_dict[key1]['comp_meth_base']) - 1].bold();
            }

        }

        if (!Object.keys(n_dict[key1]).includes('comp_meth_type')) {
            n_dict[key1]['comp_meth_type'] = "-"
        }

        var n_rec_seq = ""

        for (let i = 0; i < n_dict[key1]['rec_seq'].length; i++) {

            n_rec_seq = n_rec_seq + rec_seq_repr[i];
        }

        // console.log("Current check:", n_rec_seq, n_dict[key1]['comp_meth_type'])

        var publ_avail = [];
        var tmp_camth = camth.split("\n");

        for (let i = 0; i < tmp_camth.length; i++) {
            var tmp_publ_motif = tmp_camth[i].split("\t")[0];
            // console.log(tmp_publ_motif);
            if (n_dict[key1]['rec_seq'] == tmp_publ_motif) {
                publ_avail.push(tmp_camth[i].split("\t")[1] + " - " + tmp_camth[i].split("\t")[tmp_camth[i].split("\t").length - 1])
            }
        }

        if (n_dict[key1]['enz_type'] == "-1") {
            var ttt_enz_type = "unknown";
        } else {
            var ttt_enz_type = String(n_dict[key1]['enz_type']);
        }


        var row = table.insertRow();
        var cell0 = row.insertCell(0);
        var cell1 = row.insertCell(1);
        var cell2 = row.insertCell(2);
        var cell3 = row.insertCell(3);
        var cell4 = row.insertCell(4);
        var cell5 = row.insertCell(5);
        var cell6 = row.insertCell(6);
        var cell7 = row.insertCell(7);
        cell0.innerHTML = ttt_enz_type;
        cell1.innerHTML = n_rec_seq;
        cell2.innerHTML = n_dict[key1]['meth_base'];
        cell3.innerHTML = n_dict[key1]['meth_type'];
        cell4.innerHTML = n_dict[key1]['comp_meth_base'];
        cell5.innerHTML = n_dict[key1]['comp_meth_type'];
        cell6.innerHTML = publ_avail.join('; ');
        cell7.innerHTML = '<button type="submit" class="btn btn-primary" style="background-color: #337ab7; border-color: #2e6da4;" onclick="remove_motif(__motif_id__)">Remove</button>'.replace('__motif_id__', key1);
    }
}

function read_input_my_motifs(file_list, job_type) {
    let promises = [];
    for (let file of file_list) {
        let filePromise = new Promise(resolve => {
            let reader = new FileReader();
            reader.readAsText(file);
            reader.onload = () => resolve(reader.result);
        });
        promises.push(filePromise);
    }
    Promise.all(promises).then(fileContents => {

        var input_file = fileContents[0];


        lines = { 0: {} }
        enz_types = []
        counter = 0

        for (line of input_file.split('\n')) {
            if (line.includes('<>')) {
                counter += 1
                lines[counter] = {}
            }
            if (line.includes('enz_type')) {
                enz_types.push(Number(line.replace(/^\s+|\s+$/g, '').split('>')[1]))
            }
            if (
                (line.length > 1) &
                (!line.includes("*")) &
                (!line.includes("<>")) &
                (!line.includes("org")) &
                (!line.includes("genome")) &
                ((line.includes("rec_seq")) | (line.includes("meth")) | (line.includes("enz_type")))

            ) {
                lines[counter][line.replace(/^\s+|\s+$/g, '').split('>')[0].slice(1)] = line.replace(/^\s+|\s+$/g, '').split('>')[1]
            }
        }

        try {
            n_dict = JSON.parse($('#motifs_dict').html());
            max_key = Math.max(...Object.keys(n_dict).map(Number)) + 1;

        }
        catch (e) {
            n_dict = {};
            max_key = 1;
        }

        for (key of Object.keys(lines)) {
            if ((Object.keys(lines[key]).length != 0) & (Object.keys(lines[key]).includes('rec_seq')) & (Object.keys(lines[key]).length > 2)) {
                n_dict[Number(max_key) + Number(key)] = lines[key]
            }
        }

        console.log(n_dict);

        if ($("#rebase_placeholder").length == 0) {


            $('#motifs_dict').html(JSON.stringify(n_dict));
            $('#input_body').empty();
            $('#input_body').append(components['Submit caption']);
            $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));
            $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="upload_rebase()">Upload rebase</button>');
            $('#input_body').append('<div id="rebase_placeholder"></div>');
            $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
            draw_table_of_motifs(n_dict);
            $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="save_rebase()">Save rebase</button>');

        } else {

            $('#motifs_dict').html(JSON.stringify(n_dict));
            $('#input_body').empty();
            $('#input_body').append(components['Submit caption']);
            $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));
            $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="upload_rebase()">Upload rebase</button>');
            $('#input_body').append('<div id="rebase_placeholder"></div>');
            $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
            draw_table_of_motifs(n_dict);
            $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="save_rebase()">Save rebase</button>');
        }
    });
}

function remove_motif(key1) {

    var n_dict = JSON.parse($('#motifs_dict').html());

    delete n_dict[Number(key1)];

    $('#motifs_dict').html(JSON.stringify(n_dict));

    $('#motifs_dict').html(JSON.stringify(n_dict));
    $('#input_body').empty();
    $('#input_body').append(components['Submit caption']);
    $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));
    $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="upload_rebase()">Upload rebase</button>');
    $('#input_body').append('<div id="rebase_placeholder"></div>');
    $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
    draw_table_of_motifs(n_dict);
    $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="save_rebase()">Save rebase</button>');

}

function update_methylated_bases() {

    var curr_str = $('#tmp_input_motif_rec_seq').val();
    $('#tmp_input_motif_meth_base').empty();
    $('#tmp_input_motif_meth_base').append('<option value="-" selected>Unknown</option>');
    $('#tmp_input_motif_comp_base').empty();
    $('#tmp_input_motif_comp_base').append('<option value="-" selected>Unknown</option>');

    var rec_seq_repr = {};

    for (let i = 0; i < curr_str.length; i++) {
        rec_seq_repr[i] = curr_str.charAt(i).toLowerCase();
    }

    for (let i = 0; i < curr_str.length; i++) {

        var tt_repr = JSON.parse(JSON.stringify(rec_seq_repr));
        tt_repr[i] = tt_repr[i].toUpperCase();

        var n_rec_seq = ""

        for (let j = 0; j < curr_str.length; j++) {

            n_rec_seq = n_rec_seq + tt_repr[j];
        }
        $('#tmp_input_motif_meth_base').append(
            $('<option>', {
                value: i + 1,
                text: n_rec_seq
            }));

        $('#tmp_input_motif_comp_base').append(
            $('<option>', {
                value: i + 1,
                text: n_rec_seq
            }));

    }

}

function add_other_RM_motifs() {
    $('#rebase_placeholder').append(
        '<form>' +
        '<div class="form-row">' +
        '<div class="form-group col-md-1">' +
        '<label for="tmp_input_motif_type">Type</label>' +
        '<select class="form-select form-control" aria-label="mmeth_type" id="tmp_input_motif_type">' +
        '<option value="-1" selected>Unknown</option>' +
        '<option value="1">Type I</option>' +
        '<option value="2">Type II</option>' +
        '<option value="3">Type III</option>' +
        '<option value="4">Type IV</option>' +
        '</select>' +
        // '<input type="text" class="form-control" id="tmp_input_motif_type">' +
        '</div>' +
        '<div class="form-group col-md-3">' +
        '<label for="tmp_input_motif_rec_seq">Rec. seq.</label>' +
        '<input type="text" class="form-control" id="tmp_input_motif_rec_seq" onchange="update_methylated_bases()">' +
        '</div>' +
        '<div class="form-group col-md-2">' +
        '<label for="tmp_input_motif_meth_base">Methylated base (+)</label>' +
        '<select class="form-select form-control" aria-label="mmeth_type" id="tmp_input_motif_meth_base">' +
        '<option value="-" selected>Unknown</option>' +
        '</select>' +
        // '<input type="text" class="form-control" id="tmp_input_motif_meth_base">' +
        '</div>' +
        '<div class="form-group col-md-2">' +
        '<label for="tmp_input_motif_meth_type">Methylated base (+) type</label>' +
        '<select class="form-select form-control" aria-label="mmeth_type" id="tmp_input_motif_meth_type">' +
        '<option value="-" selected>Unknown</option>' +
        '<option value="6">m6A</option>' +
        '<option value="4">m4C</option>' +
        '<option value="5">m5C</option>' +
        '</select>' +
        // '<input type="text" class="form-control" id="tmp_input_motif_meth_type">' +
        '</div>' +
        '<div class="form-group col-md-2">' +
        '<label for="tmp_input_motif_comp_base">Methylated base (-)</label>' +
        '<select class="form-select form-control" aria-label="mmeth_type" id="tmp_input_motif_comp_base">' +
        '<option value="-" selected>Unknown</option>' +
        '</select>' +
        // '<input type="text" class="form-control" id="tmp_input_motif_comp_base">' +
        '</div>' +
        '<div class="form-group col-md-2">' +
        '<label for="tmp_input_motif_comp_type">Methylated base (-) type</label>' +
        '<select class="form-select form-control" aria-label="mmeth_type" id="tmp_input_motif_comp_type">' +
        '<option value="-" selected>Unknown</option>' +
        '<option value="6">m6A</option>' +
        '<option value="4">m4C</option>' +
        '<option value="5">m5C</option>' +
        '</select>' +
        // '<input type="text" class="form-control" id="tmp_input_motif_comp_type">' +
        '</div>' +
        '</div>' +
        '<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4; margin-bottom:0.5%" id="add_other_RM_motifs_to_the_table_button" onclick="add_other_RM_motifs_to_the_table()">Add motif</button>' +
        '</form>'
    )
    // $('#rebase_placeholder').append();
    $('#__id_add_icon__').remove();
}

function add_other_RM_motifs_to_the_table() {
    var tmp_input_motif_type = $('#tmp_input_motif_type').val();
    console.log("tmp_input_motif_type", tmp_input_motif_type);
    // console.log(tmp_input_motif_type == "");
    var tmp_input_motif_rec_seq = $('#tmp_input_motif_rec_seq').val();
    console.log("tmp_input_motif_rec_seq", tmp_input_motif_rec_seq);
    // console.log(tmp_input_motif_rec_seq == "");
    var tmp_input_motif_meth_base = $('#tmp_input_motif_meth_base').val();
    console.log("tmp_input_motif_meth_base", tmp_input_motif_meth_base);
    // console.log(tmp_input_motif_meth_base == "");
    var tmp_input_motif_meth_type = $('#tmp_input_motif_meth_type').val();
    console.log("tmp_input_motif_meth_type", tmp_input_motif_meth_type);
    // console.log(tmp_input_motif_meth_type == "");
    var tmp_input_motif_comp_base = $('#tmp_input_motif_comp_base').val();
    console.log("tmp_input_motif_comp_base", tmp_input_motif_comp_base);
    // console.log(tmp_input_motif_comp_base == "");
    var tmp_input_motif_comp_type = $('#tmp_input_motif_comp_type').val();
    console.log("tmp_input_motif_comp_type", tmp_input_motif_comp_type);
    // console.log(tmp_input_motif_comp_type == "");

    var letters = [...new Set(tmp_input_motif_rec_seq)]
    letters.sort()
    var letters = letters.join('').toUpperCase();
    n1 = "ABCDGHKMNRSTVWY".split('');
    n2 = letters.split('');

    var correct = true;

    if (tmp_input_motif_rec_seq == "") {
        correct = false;
        alert("Input rec. sec. not defined!")
    }

    if (!n2.every(i => n1.includes(i))) {
        correct = false;
        alert("Rec. seq. not valid!");
    }

    if (tmp_input_motif_type != "") {
        if (!["-1", "1", "2", "3", "4"].includes(String(tmp_input_motif_type))) {
            correct = false;
            alert("Type not valid!");
        }
    }

    if (tmp_input_motif_meth_base != "-") {
        if (!(((Number(tmp_input_motif_meth_base) <= tmp_input_motif_rec_seq.length) & (Number(tmp_input_motif_meth_base) > 0)))) {
            correct = false;
            alert("Invalid methylated base!");
        }
    }

    if (((tmp_input_motif_meth_type != '-') & (tmp_input_motif_meth_type != ""))) {
        if (!(["-", "4", "5", "6"].includes(String(tmp_input_motif_meth_type)))) {
            correct = false;
            alert("Invalid methylated base type!");
        }
    }

    if (tmp_input_motif_comp_base != "-") {
        if (!(((Number(tmp_input_motif_comp_base) <= tmp_input_motif_rec_seq.length) & (Number(tmp_input_motif_comp_base) > 0)))) {
            // if (!(Number(tmp_input_motif_comp_base) <= tmp_input_motif_rec_seq.length)) {
            correct = false;
            alert("Invalid methylated base!");
        }
    }

    if (tmp_input_motif_comp_type != "") {
        if (!(["-", "4", "5", "6"].includes(String(tmp_input_motif_comp_type)))) {
            correct = false;
            alert("Invalid methylated base type!");
        }
    }
    if (correct) {

        var n_dict = JSON.parse($('#motifs_dict').html());

        if (Object.keys(n_dict).length > 0) {
            var obj_keys = Math.max(...Object.keys(n_dict));
        } else {
            var obj_keys = 0;
        }

        if (tmp_input_motif_type == "") {
            var tmp_input_motif_type = '-1'
        }
        if (tmp_input_motif_meth_base == "") {
            var tmp_input_motif_meth_base = '-'
        }
        if (tmp_input_motif_meth_type == "") {
            var tmp_input_motif_meth_type = '-'
        }
        if (tmp_input_motif_comp_base == "") {
            var tmp_input_motif_comp_base = '-'
        }
        if (tmp_input_motif_comp_type == "") {
            var tmp_input_motif_comp_type = '-'
        }

        n_dict[obj_keys + 1] = {
            "enz_type": tmp_input_motif_type,
            "rec_seq": tmp_input_motif_rec_seq,
            "meth_base": tmp_input_motif_meth_base,
            "meth_type": tmp_input_motif_meth_type,
            "comp_meth_base": tmp_input_motif_comp_base,
            "comp_meth_type": tmp_input_motif_comp_type
        };// alert(JSON.stringify(n_dict));

        // $("#rebase_output_table").remove();

        // if ($("#rebase_output_table").length > 0) {
        //     $("#rebase_output_table").remove();
        //     $('#upload_rebase_button').remove();
        // }

        // draw_table_of_motifs(n_dict);

        // alert("Here1");

        $('#motifs_dict').html(JSON.stringify(n_dict));

        $('#motifs_dict').html(JSON.stringify(n_dict));
        $('#input_body').empty();
        $('#input_body').append(components['Submit caption']);
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));
        $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="upload_rebase()">Upload rebase</button>');
        $('#input_body').append('<div id="rebase_placeholder"></div>');
        $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
        draw_table_of_motifs(n_dict);
        $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="save_rebase()">Save rebase</button>');

        // if ($("#rebase_placeholder").length > 0) {

        //     alert("Here4");

        //     $('#rebase_placeholder').remove();
        //     $('#input_body').append('<div id="rebase_placeholder"></div>');
        //     $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
        //     $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="save_rebase_button" onclick="save_rebase()">Save rebase</button>');

        // }

        // $('#upload_rebase_button').remove();
        // $('#rebase_placeholder').remove();
        // $('#upload_input_file').remove();
        // $('#mymotif_choice').remove();


        // $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="save_rebase()">Save rebase</button>');     
    } else {

        var n_dict = JSON.parse($('#motifs_dict').html());

        // if ($("#rebase_output_table").length > 0) {
        //     $("#rebase_output_table").remove();
        //     $('#upload_rebase_button').remove();
        // }

        // draw_table_of_motifs(n_dict);

        // alert("Here2");

        $('#motifs_dict').html(JSON.stringify(n_dict));

        $('#motifs_dict').html(JSON.stringify(n_dict));
        $('#input_body').empty();
        $('#input_body').append(components['Submit caption']);
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));
        $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="upload_rebase()">Upload rebase</button>');
        $('#input_body').append('<div id="rebase_placeholder"></div>');
        $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
        draw_table_of_motifs(n_dict);
        $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="save_rebase()">Save rebase</button>');

        // if ($("#rebase_placeholder").length > 0) {

        //     alert("Here4");

        //     $('#rebase_placeholder').remove();
        //     $('#input_body').append('<div id="rebase_placeholder"></div>');
        //     $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
        //     $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="save_rebase_button" onclick="save_rebase()">Save rebase</button>');

        // }

        // $("#rebase_output_table").remove();

        // draw_table_of_motifs(n_dict);

        // $('#motifs_dict').html(JSON.stringify(n_dict));

        // $('#upload_rebase_button').remove();
        // $('#rebase_placeholder').remove();
        // // $('#upload_input_file').remove();
        // // $('#mymotif_choice').remove();
        // $('#input_body').append('<div id="rebase_placeholder"></div>');
        // $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM motif').replace('__action__', 'add_other_RM_motifs()'));
        // // $('#__id_add_icon__').show();
        // $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="save_rebase()">Save rebase</button>');     
    }

}

function save_rebase() {
    var n_dict = JSON.parse($('#motifs_dict').html()); //alert($('#motifs_dict').html());

    var properties = { type: 'text/plain' };

    var out_file = [];

    var type_iv_present = false;

    console.log("Saving rebase:", n_dict)

    for (key of Object.keys(n_dict)) {

        if (n_dict[key]["enz_type"] != undefined) {
            out_file.push("<enz_type>" + n_dict[key]["enz_type"] + "\n");

            if (n_dict[key]["enz_type"] == "4") {
                // alert(n_dict[key]["enz_type"]);
                type_iv_present = true;
            }

        }

        out_file.push("<rec_seq>" + n_dict[key]["rec_seq"] + "\n");

        if (Object.keys(n_dict[key]).includes("meth_base")) {
            out_file.push("<meth_base>" + n_dict[key]["meth_base"] + "\n");
        }

        if (Object.keys(n_dict[key]).includes("meth_type")) {
            out_file.push("<meth_type>" + n_dict[key]["meth_type"] + "\n");
        }

        if (Object.keys(n_dict[key]).includes("comp_meth_base")) {
            out_file.push("<comp_meth_base>" + n_dict[key]["comp_meth_base"] + "\n");
        }

        if (Object.keys(n_dict[key]).includes("comp_meth_type")) {
            out_file.push("<comp_meth_type>" + n_dict[key]["comp_meth_type"] + "\n");
        }
        out_file.push("<>\n");
    }

    if (type_iv_present) {
        if (!confirm('Type IV RM present in the target strain, do you want to confirm?')) {
            return
        }
    } else {
        if (!confirm('Type IV RM not present in the target strain, do you want to confirm?')) {
            return
        }
    }

    try {
        var file = new File(out_file, "file.txt", properties);
    } catch (e) {
        var file = new Blob(out_file, properties);
    }

    download(file, "sytogen_suite_my_motifs_out.txt")

}

function upload_rebase() {

    var input_file = [$('#output_rebase')[0].files[0]];

    read_input_my_motifs(input_file);
}