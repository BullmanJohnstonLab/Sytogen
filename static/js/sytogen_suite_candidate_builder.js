function trigger_current_slice(key) {

    var tmp_current_considered_slice = $("#current_considered_slice").value;

    if ((tmp_current_considered_slice != null) ) {


        if ($('#' + tmp_current_considered_slice).getAttribute("selected") != "+1") {
            $('#' + tmp_current_considered_slice).setAttribute('selected', "-1");
            $('#' + tmp_current_considered_slice).style.color = "black";
            $('#' + tmp_current_considered_slice).style.background = "white";

            if ( $("#configFile").value != null) {
                var tmp_config_file = JSON.parse( $("#configFile").value );
                delete tmp_config_file[$('#' + tmp_current_considered_slice).value];
                $("#configFile").value = JSON.stringify(tmp_config_file);
                $("#current_considered_slice").value = key;
            }
        }
    }

    if ( document.getElementById(key).getAttribute('selected') == "-1" ) {
        document.getElementById(key).setAttribute('selected', "0");
        document.getElementById(key).style.color = "white";
        document.getElementById(key).style.background = "#d8b365";

        document.getElementById('current_considered_slice').value = key;
        document.getElementById('current_alternative_slice').value = null;
        retrieve_list_of_candidate_alternatives(key);
        
    } else if (document.getElementById(key).getAttribute('selected') == "0") {
        document.getElementById(key).setAttribute('selected', "-1");
        document.getElementById(key).style.color = "black";
        document.getElementById(key).style.background = "white";

        if (document.getElementById("configFile").value != null) {
            var tmp_config_file = JSON.parse(document.getElementById("configFile").value);
            delete tmp_config_file[key];
            document.getElementById("configFile").value = JSON.stringify(tmp_config_file);
            document.getElementById('current_considered_slice').value = null;
        }

        document.getElementById("temporary_select_slice").innerHTML = null;
        document.getElementById('current_alternative_slice').value = null;

    } else if (document.getElementById(key).getAttribute('selected') == "+1") {
        
        if (document.getElementById('current_alternative_slice').value != null) {
            document.getElementById('current_alternative_slice').getAttribute("selected", "-1");
            document.getElementById('current_alternative_slice').style.color = "black";
            document.getElementById('current_alternative_slice').style.background = "white";
        }

        // console.log(document.getElementById("configFile").value);

        document.getElementById(key).setAttribute('selected', "+1");
        document.getElementById(key).style.color = "white";
        document.getElementById(key).style.background = "#5ab4ac";

        var tmp_config_file = JSON.parse(document.getElementById("configFile").value);
        document.getElementById('current_alternative_slice').value = tmp_config_file[key];

        retrieve_list_of_candidate_alternatives(key);
        
    }
    
    return false;
}

function trigger_current_candidate_alternative(key) {

    // console.log(document.getElementById('current_considered_slice').value);
    // console.log(document.getElementById('current_alternative_slice').value);
    // console.log(document.getElementById(key).getAttribute('selected'));

    if (document.getElementById(key).getAttribute('selected') == "-1") {


        document.getElementById('current_alternative_slice').value = key;
        retrieve_list_of_candidate_alternatives(document.getElementById('current_considered_slice').value);
        document.getElementById(key).setAttribute('selected', "0");
        // document.getElementById(key).style.color = "white";
        // document.getElementById(key).style.background = "#5ab4ac";

        document.getElementById(document.getElementById('current_considered_slice').value).setAttribute('selected', "+1");
        document.getElementById(document.getElementById('current_considered_slice').value).style.color = "white";
        document.getElementById(document.getElementById('current_considered_slice').value).style.background = "#5ab4ac";

        if (document.getElementById("configFile").value != null) {
            var tmp_config_file = JSON.parse(document.getElementById("configFile").value);
            tmp_config_file[document.getElementById('current_considered_slice').value] = key;
            document.getElementById("configFile").value = JSON.stringify(tmp_config_file);
        } else {
            var tmp_config_file = {};
            tmp_config_file[document.getElementById('current_considered_slice').value] = key;
            document.getElementById("configFile").value = JSON.stringify(tmp_config_file);
        }

    } else if (document.getElementById(key).getAttribute('selected') == "0") {
        
        // console.log("HERE")

        document.getElementById(key).setAttribute('selected', "-1");
        document.getElementById('current_alternative_slice').value = null;
        retrieve_list_of_candidate_alternatives(document.getElementById('current_considered_slice').value);
        // document.getElementById(key).style.color = "black";
        // document.getElementById(key).style.background = "white";
        // console.log("I AM HERE!")
        document.getElementById(document.getElementById('current_considered_slice').value).setAttribute('selected', "0");
        document.getElementById(document.getElementById('current_considered_slice').value).style.color = "white";
        document.getElementById(document.getElementById('current_considered_slice').value).style.background = "#d8b365";

        if (document.getElementById("configFile").value != null) {
            var tmp_config_file = JSON.parse(document.getElementById("configFile").value);
            tmp_config_file[document.getElementById('current_considered_slice').value] = null;
            document.getElementById("configFile").value = JSON.stringify(tmp_config_file);
        } else {
            var tmp_config_file = {};
            tmp_config_file[document.getElementById('current_considered_slice').value] = key;
            document.getElementById("configFile").value = JSON.stringify(tmp_config_file);
        }
    }


    // console.log(document.getElementById('current_considered_slice').value);
    // console.log(document.getElementById('current_alternative_slice').value);
    // console.log(document.getElementById(key).getAttribute('selected'));

    // console.log(document.getElementById("configFile").value);
    return false;
};

function retrieve_list_of_slices() {

    $("#retrieve_list_of_slices_button").remove();
    $('#input_body').append(components['Loader']);
    
    var data = new FormData();
    var job_id = $("#syn_cand_builder_job_id").val();
    data.append("job_id", job_id);

    $('#current_job_id').val(job_id);

    fetch(proxy_prefix + '/get_synpl_list/', {
        method: 'POST',
        body: data
    }).then(function (response) {
        return response.text();
    }).then(function (text) {
        if (text.includes(".synpl")) {

            tmp_var_text = JSON.parse(text);
            var keys = Object.keys(tmp_var_text);
            var list_slices_options = document.getElementById("select_slice");
            
            document.getElementById('synpl_slice_list').value = keys;
            for (i = 0; i < keys.length; i++) {
                var a = document.createElement("a");
                a.appendChild(document.createElement('br'));
                a.appendChild(document.createTextNode(keys[i]));
                for (j = 0; j < tmp_var_text[keys[i]].length; j++) {
                    a.appendChild(document.createElement('br'));
                    a.appendChild(document.createTextNode(tmp_var_text[keys[i]][j].split(" ").join(".")));
                }
                a.setAttribute("class", "list-group-item list-group-item-action");
                a.setAttribute("id", keys[i]);
                a.setAttribute("selected", "-1");
                a.addEventListener("click", trigger_current_slice.bind(null, keys[i]));
                
                list_slices_options.appendChild(a);
            }
            $('#list_of_synpl_files_box').css('visibility', 'visible');
            $('#loader').remove();
            $('#input_body').append(components['Candidate syngenic oth input']);
            $('#input_body').append(components['Candidate syngenic oth button']);
        } else {
            alert("No available output!");
        }
    });
    return false;
}


function retrieve_list_of_candidate_alternatives(text) {

    $('#input_body').append(components['Loader']);
    $('#loader').css('float', 'right');
    // $('#build_candidate_syngenic_button_').remove();    

    $('#list_of_alternative_files_box').attr('style', 'visibility: visible');
    $('#other_alternatives').attr('style', 'visibility: visible');
    $('#build_button_box').attr('style', 'visibility: visible;');


    document.getElementById("current_considered_slice").value = text;

    var data = new FormData();
    var job_id = $("#syn_cand_builder_job_id").val();
    data.append("job_id", job_id);
    data.append("slice_synpl", text)
    fetch(proxy_prefix + '/single_candidate_selection/', {
        method: 'POST',
        body: data
    }).then(function (response) {
        return response.text();
    }).then(function (text) {

        tmp_var_text = JSON.parse(text);

        var keys = Object.keys(tmp_var_text);

        var list_slices_options = document.getElementById("temporary_select_slice");
        while (list_slices_options.firstChild)
            list_slices_options.removeChild(list_slices_options.firstChild);

        for (i = 0; i < keys.length; i++) {
            var a = document.createElement("a");
            a.appendChild(document.createElement('br'));
            a.appendChild(document.createTextNode(keys[i]));
            for (j = 0; j < tmp_var_text[keys[i]].length; j++) {
                a.appendChild(document.createElement('br'));
                a.appendChild(document.createTextNode(tmp_var_text[keys[i]][j].split(" ").join(".")));
            }
            a.setAttribute("class", "list-group-item list-group-item-action");
            a.setAttribute("id", keys[i]);
            a.setAttribute("selected", "-1");
            a.addEventListener("click", trigger_current_candidate_alternative.bind(null, keys[i]));
            list_slices_options.appendChild(a);
        }

        if (document.getElementById('current_alternative_slice').value != null) {
            document.getElementById(document.getElementById('current_alternative_slice').value).setAttribute("selected", "0");
            document.getElementById(document.getElementById('current_alternative_slice').value).style.color = "white";
            document.getElementById(document.getElementById('current_alternative_slice').value).style.backgroundColor = "#5ab4ac";
        }
        // $('#input_body').append(components['Candidate syngenic oth button']);
        $('#loader').remove();
    });
    return "END";
};

function build_candidate_syngenic() {

    

    var data = new FormData();
    var job_id = $("#syn_cand_builder_job_id").val();
    var configFile = document.getElementById('configFile').value;

    var all_keys = document.getElementById('synpl_slice_list').value;

    console.log(configFile);

    if (configFile != undefined) {

        var tmp_json = JSON.parse(configFile);

        console.log(Object.values(tmp_json));
        console.log([... new Set(Object.values(tmp_json))]);
        console.log(new Set(Object.values(tmp_json)))

        console.log(String([... new Set(Object.values(tmp_json))]) == String([null]))
    }

    if (configFile == undefined) {
        alert('No candidate alternative slice selected for any of the generated slices!');
    } else if (configFile != undefined) {
        var tmp_configFile = JSON.parse(configFile);
        if (String([... new Set(Object.values(tmp_configFile))]) == String([null])) {
            alert('No candidate alternative slice selected for any of the generated slices!')
        } else {

            $('#input_body').append(components['Loader']);
            $('#loader').css('float', 'right');
            $('#build_candidate_syngenic_button_').remove();

            for (i = 0; i < all_keys.length; i++) {
                if (!(tmp_configFile.hasOwnProperty(all_keys[i]))) {
                    tmp_configFile[all_keys[i]] = document.getElementById('other_alternatives_input').value;
                }
            }

            var configFile = JSON.stringify(tmp_configFile);

            data.append("configFile", configFile);
            data.append("job_id", job_id);

            fetch(proxy_prefix + '/build_candidate_syngenic/', {
                method: 'POST',
                body: data
            }).then((data) => {

                $("#build_candidate_syngenic_button_").remove();

                $('#loader').remove();

                // $('#input_body').append('Your job id is: '.bold() + $('#current_job_id').val());

                $('#input_body').append(components['Check status button']);

                $('#current_job_type').val('candidate_syngenic');

            });
        } 
    }
}