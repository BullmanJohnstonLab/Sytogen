function load_page_components(current_page) {
    if ((current_page == 'home_page') | (current_page == null)) {
        $("#current_central_panel").load(proxy_prefix + "/static/documentation/Welcome.html");
    } else if (current_page == 'data_samples') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Data samples'));
        $('#current_central_panel').append("<h4>RM systems information (config file from REBASE)</h4>");
        $('#current_central_panel').append("<p class=\"data_sample\">" + rm_target_motifs_s_aureus + "</p>");
        $('#current_central_panel').append("<h4>Input genetic tool (GenBank format)</h4>");
        $('#current_central_panel').append("<p class=\"data_sample\">" + curr_pl + "</p>");
        $('#current_central_panel').append("<h4>Input genetic tool backbone (Fasta format)</h4>");
        $('#current_central_panel').append("<p class=\"data_sample\">" + curr_backbone + "</p>");
    } else if (current_page == 'documentation') {
        // $('#current_central_panel').append("<zero-md src='"+proxy_prefix+"/static/documentation/README.md'></zero-md>")
        $("#current_central_panel").load(proxy_prefix + "/static/documentation/Documentation.html");
    } else if (current_page.includes('take_a_tour')) {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'SyToGen Suite - Tour'));
        $('#current_central_panel').append(components['Page subtitle'].replace('__VAR_STR__', 'The Sytogen suite'));
        if (current_page == 'take_a_tour') {            
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Start the tour']));
            $('#take_a_tour_sytogen_suite').addClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('sytogen_suite')) {
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Start the tour']));
            $('#take_a_tour_sytogen_suite').addClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('codon_bias_estimator')) {            
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Tour codon bias estimator']));
            $('#take_a_tour_sytogen_suite').removeClass('active');
            $('#take_a_tour_codon_bias_estimator').addClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('synplest')) {            
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Tour synplest']));
            $('#take_a_tour_sytogen_suite').removeClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').addClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('candidate_syngenic')) {            
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Tour candidate syngenic']));
            $('#take_a_tour_sytogen_suite').removeClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').addClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('check_status')) {
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Tour check status']));
            $('#take_a_tour_sytogen_suite').removeClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').addClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('check_input')) {            
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Tour check input']));
            $('#take_a_tour_sytogen_suite').removeClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').addClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('check_output')) {            
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Tour check output']));
            $('#take_a_tour_sytogen_suite').removeClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').addClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('tool_representation')) {            
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Tour tool representation']));
            $('#take_a_tour_sytogen_suite').removeClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').addClass('active');
            $('#take_a_tour_tool_preprocess').removeClass('active');
        } else if (current_page.includes('tool_preprocess')) {            
            $('#current_central_panel').append(components['Tour card'].replace('__BODY__', components['Tour tool preprocess']));
            $('#take_a_tour_sytogen_suite').removeClass('active');
            $('#take_a_tour_codon_bias_estimator').removeClass('active');
            $('#take_a_tour_synplest').removeClass('active');
            $('#take_a_tour_candidate_syngenic').removeClass('active');
            $('#take_a_tour_check_status').removeClass('active');
            $('#take_a_tour_check_input').removeClass('active');
            $('#take_a_tour_check_output').removeClass('active');
            $('#take_a_tour_tool_representation').removeClass('active');
            $('#take_a_tour_tool_preprocess').addClass('active');
        }
    } else if (current_page == 'contacts') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'SyToGen Suite - Contacts'));
        $('#current_central_panel').append(components['contacts']);
    } else if (current_page == 'mymotifs') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'MyMotifs'));
        $('#current_central_panel').append(components['Tour MyMotifs']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append(components['Submit caption']);

        $('#input_body').append(components['MyMotif choice']);

        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));

        $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="upload_rebase_button" onclick="upload_rebase()">Upload rebase</button>');

        
    } else if (current_page == 'motifsfinder') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'MotifsFinder'));
        $('#current_central_panel').append(components['Tour MotifsFinder']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append('<div id="bug_placeholder""></div>');
        $('#bug_placeholder').append(components['Submit caption']);
        $('#bug_placeholder').append("<p>Input genetic tool(s) must be in GenBank format.</p>");
        $('#bug_placeholder').append(components["Add_icon"].replace('__info__', 'Add other genetic tools').replace('__action__', 'add_other_bugs()'));
        // $('#bug_placeholder').append(components["Bug_info"].replaceAll('__bug_genus__', '1__input_genetic_tool_genus').replaceAll('__bug_species__', '1__input_genetic_tool_species').replaceAll('__bug_strain__', '1__input_genetic_tool_strain'));
        $('#bug_placeholder').append(components['Input file'].replace('__VAR_LABEL__', 'Input gen. tool').replace('__VAR_ID__', '1__input_genetic_tool_file'));
        //
        $('#input_body').append('<hr>')
        // 
        // bug_ids, rebase_ids
        $('#input_body').append('<div id="rebase_placeholder"></div>');
        $('#rebase_placeholder').append(components["Add_icon"].replace('__info__', 'Add other RM systems info (REBASE details)').replace('__action__', 'add_other_rebase()'));
        // $('#rebase_placeholder').append(components["Bug_info"].replaceAll('__bug_genus__', '1__input_rebase_genus').replaceAll('__bug_species__', '1__input_rebase_species').replaceAll('__bug_strain__', '1__input_rebase_strain'));
        $('#rebase_placeholder').append(components['Input file'].replace('__VAR_LABEL__', 'Input REBASE').replace('__VAR_ID__', '1__input_rebase_file'));

        $('#bug_ids').html(JSON.stringify(['1']));
        $('#rebase_ids').html(JSON.stringify(['1']));

        //
        $('#input_body').append(components['Submit button']);

    } else if (current_page == 'codon_bias_estimator') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Codon bias estimator'));
        $('#current_central_panel').append(components['Tour codon bias estimator']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append(components['Submit caption']);
        $('#input_body').append("<p>Input genome must be in GenBank format.</p>");
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'Input genome').replace('__VAR_ID__', 'input_genome'));
        $('#input_body').append(components['Codon table choice']);
        $('#input_body').append(components['Submit button']);
    } else if (current_page == 'synplest') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Syngeneic Tool Estimator'));
        $('#current_central_panel').append(components['Tour synplest']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append(components['Submit caption']);
        $('#input_body').append("<p>Input tool and genome must be in GenBank format.</p>");
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'Input plasmid').replace('__VAR_ID__', 'input_plasmid'));
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'Input genome').replace('__VAR_ID__', 'input_genome'));
        $('#input_body').append(components['Codon table choice']);
        $('#input_body').append(components['Remove backbone choice']);
        $('#input_body').append(components['Submit button']);
    } else if (current_page == 'check_status') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Check job status'));
        $('#current_central_panel').append(components['Tour check status']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append(components['Input job id'].replace('__JOB_ID__', 'job_id_check'));
        $('#input_body').append(components['Check status button']);
    } else if (current_page == 'candidate_syngenic') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Candidate syngenic builder'));
        $('#current_central_panel').append(components['Tour candidate syngenic']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append(components['Candidate syngenic input']);
        $('#input_body').append(components['Retrieve generated slices button']);
        $('#current_central_panel').append(components['Candidate syngenic']);
    } else if (current_page == 'check_input') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Check input'));
        $('#current_central_panel').append(components['Tour check input']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append("<p>Input sequence must be in GenBank format.</p>");
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'Input sequence').replace('__VAR_ID__', 'input_sequence'));
        $('#input_body').append(components['Submit button']);
    } else if (current_page == 'check_output') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Check output'));
        $('#current_central_panel').append(components['Tour check output']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append(components['Submit caption']);
        $('#input_body').append("<p>Input tools must be in GenBank format.</p>");
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'Original tool').replace('__VAR_ID__', 'or_input_plasmid'));
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'Syngenic tool').replace('__VAR_ID__', 'nw_input_plasmid'));
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));
        $('#input_body').append(components['Codon table choice']);
        $('#input_body').append(components['Submit button']);
    } else if (current_page == 'tool_representation') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Tool representation'));
        $('#current_central_panel').append(components['Tour tool representation']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append(components['Submit caption']);
        $('#input_body').append("<p>Input tool must be in GenBank format.</p>");
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'Input tool').replace('__VAR_ID__', 'input_sequence'));
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'RM information').replace('__VAR_ID__', 'output_rebase'));
        $('#input_body').append(components['Submit button']);
    } else if (current_page == 'tool_preprocess') {
        $('#current_central_panel').append(components['Page title'].replace('__VAR_STR__', 'Tool preprocess'));
        $('#current_central_panel').append(components['Tour tool preprocess']);
        $('#current_central_panel').append(components['Input box']);
        $('#input_body').append(components['Submit caption']);
        $('#input_body').append("<p>Input tool must be in GenBank format.</p>");
        $('#input_body').append(components['Input file'].replace('__VAR_LABEL__', 'Input tool').replace('__VAR_ID__', 'input_sequence'));
        $('#input_body').append(components['Caption'].replace('__VAR_CAPTION__', 'If you decide to remove the backbone, it will not be considered by the Syngeneic Tool Estimator. <br> WARNING: When the backbone is removed, the new proposed sequence will be only the (circular) sequence of the mini-circle. <br>In this case you will need to add again the backbone sequence, if required.').replace('__CAP_ID__', 'Remove_backbone_caption'));
        $('#input_body').append("<p>Backbone sequence must be in Fasta format.</p>");
        $('#input_body').append(components['Backbone sequence'].replace('__BACKBONE_SEQ__', 'backbone_seq_box'));
        $('#input_body').append(components['Submit button']);
    }
}

function load_current_page() {
    var path = window.location.href;

    var current_page = path.split("#")[1];

    $("#layoutSidenav_nav .sb-sidenav a.nav-link").each(function () {
        if (this.href === path) {
            $(this).addClass("active");
        } else {
            $(this).removeClass("active");
        }
    });

    $('#current_central_panel').empty();
    $('#current_considered_slice').empty();
    $('#current_alternative_slice').empty();
    $('#synpl_slice_list').empty();
    $('#configFile').empty();
    $('#selection_object').empty();
    $('#current_job_id').empty();
    $('#current_job_type').empty();
    load_page_components(current_page);

}

function add_other_bugs() {
    
    var curr_ids = JSON.parse($('#bug_ids').html());
    var last_id = Number(curr_ids[curr_ids.length - 1]) + 1;
    curr_ids.push(last_id);
    $('#bug_ids').html(JSON.stringify(curr_ids));

    // $('#bug_placeholder').append(components["Bug_info"].replaceAll('__bug_genus__', String(last_id) + '__input_genetic_tool_genus').replaceAll('__bug_species__', String(last_id) + '__input_genetic_tool_species').replaceAll('__bug_strain__', last_id + '__input_genetic_tool_strain'));
    $('#bug_placeholder').append(components['Input file'].replace('__VAR_LABEL__', 'Input gen. tool').replace('__VAR_ID__', String(last_id) + '__input_genetic_tool_file'));

}

function add_other_rebase() {
    
    var curr_ids = JSON.parse($('#rebase_ids').html());
    var last_id = Number(curr_ids[curr_ids.length - 1]) + 1;
    curr_ids.push(last_id);
    $('#rebase_ids').html(JSON.stringify(curr_ids));

    // $('#rebase_placeholder').append(components["Bug_info"].replaceAll('__bug_genus__', String(last_id) + '__input_rebase_genus').replaceAll('__bug_species__', String(last_id) + '__input_rebase_species').replaceAll('__bug_strain__', String(last_id) + '__input_rebase_strain'));
    $('#rebase_placeholder').append(components['Input file'].replace('__VAR_LABEL__', 'Input REBASE').replace('__VAR_ID__', String(last_id) + '__input_rebase_file'));

}



function get_backbone_choice(this_obj) {

    var choice = $('input[name="optradio"]:checked').val()

    if (choice == 'Consider the complete sequence') {

        $('#Remove_backbone_caption').remove();
        $('#backbone_seq_box').remove();

        $('#submit_button').remove();
        $('#input_body').append(components['Submit button']);

        if ($('#remove_backbone').length) {
            $('#remove_backbone').html('false');
        } else {
            $('#main_body').append('<div id="remove_backbone" style="display: none"></div>');
            $('#remove_backbone').html('false');

        }

    } else if (choice == 'Remove backbone') {

        $('#get_target_position_choice').remove();
        $('#submit_button').remove();
        $('#input_body').append(components['Caption'].replace('__VAR_CAPTION__', 'If you decide to remove the backbone, it will not be considered by the Syngeneic Tool Estimator. <br> WARNING: When the backbone is removed, the new proposed sequence will be only the (circular) sequence of the mini-circle. <br>In this case you will need to add again the backbone sequence, if required.').replace('__CAP_ID__', 'Remove_backbone_caption'));
        $('#input_body').append("<p>Backbone sequence must be in Fasta format.</p>");
        $('#input_body').append(components['Backbone sequence'].replace('__BACKBONE_SEQ__', 'backbone_seq_box'));
        $('#input_body').append(components['Submit button']);

        if ($('#remove_backbone').length) {
            $('#remove_backbone').html('true');
        } else {
            $('#main_body').append('<div id="remove_backbone" style="display: none"></div>');
            $('#remove_backbone').html('true');
        }

    }
}

function submit_job() {

    $("#submit_button").remove();

    $('#input_body').append(components['Loader']);

    var path = window.location.href;
    var current_page = path.split("#")[1];
    
    if (current_page == 'codon_bias_estimator') {

        var fd = new FormData();
        var input_genome = $('#input_genome')[0].files[0];
        var codon_table = $("#codon_table").val();

        fd.append('input_genome', input_genome);
        fd.append('codon_table', codon_table);

        $.ajax({
            type: "POST",
            url: proxy_prefix + "/run_codon_bias/",
            data: fd,
            processData: false,
            contentType: false,
            dataType: "json",
            success: function (data) {

                if (data['type'] == 'Success') {
                    $('#output_status').remove();

                    $('#loader').remove();

                    $('#input_body').append('Your job id is: '.bold() + data.job_id);
                    
                    $('#input_body').append(components['Check status button']);

                    $('#current_job_id').val(data.job_id);
                    $('#current_job_type').val('codon_bias');

                } else {
                    alert(data.message);
                    $('#loader').remove();
                }
            },
            error: function (req, status, err) {
                var resp = req.responseJSON;
            }
        });

    } else if (current_page == 'synplest') {

        var fd = new FormData();

        var input_plasmid = $('#input_plasmid')[0].files[0];
        var output_rebase = $('#output_rebase')[0].files[0];
        var input_genome = $('#input_genome')[0].files[0];
        var codon_table = $("#codon_table").val();
        var remove_backbone = $('#remove_backbone').html();
        var backbone_input_sequence = $('#backbone_input_sequence').val();

        fd.append('input_plasmid', input_plasmid);
        fd.append('output_rebase', output_rebase);
        fd.append('input_genome', input_genome);
        fd.append('codon_table', codon_table);
        fd.append('remove_backbone', remove_backbone);
        fd.append('backbone_input_sequence', backbone_input_sequence);

        $.ajax({
            type: "POST",
            url: proxy_prefix + "/run_synplest/",
            data: fd,
            processData: false,
            contentType: false,
            dataType: "json",
            success: function (data) {

                if (data['type'] == 'Success') {
                    $('#output_status').remove();

                    $('#loader').remove();

                    $('#input_body').append('Your job id is: '.bold() + data.job_id);

                    $('#input_body').append(components['Check status button']);

                    $('#current_job_id').val(data.job_id);
                    $('#current_job_type').val('synplest');

                    $('input[name="optradio"]').prop("disabled", true);

                } else {
                    alert(data.message);
                    $('#loader').remove();
                }

            },
            error: function (req, status, err) {
                var resp = req.responseJSON;
            }
        });

    } else if (current_page == 'check_input') {
        var fd = new FormData();

        var input_sequence = $('#input_sequence')[0].files[0];

        fd.append('input_sequence', input_sequence);

        $.ajax({
            type: "POST",
            url: proxy_prefix + "/run_check_input/",
            data: fd,
            processData: false,
            contentType: false,
            dataType: "json",
            success: function (data) {

                if (data['type'] == 'Success') {
                    $('#output_status').remove();

                    $('#loader').remove();

                    $('#input_body').append('Your job id is: '.bold() + data.job_id);

                    $('#input_body').append(components['Check status button']);

                    $('#current_job_id').val(data.job_id);
                    $('#current_job_type').val('check_input');

                } else {
                    alert(data.message);
                    $('#loader').remove();
                }
            },
            error: function (req, status, err) {
                var resp = req.responseJSON;
            }
        });
    } else if (current_page == 'check_output') {
        var fd = new FormData();

        var or_input_plasmid = $('#or_input_plasmid')[0].files[0];
        var nw_input_plasmid = $('#nw_input_plasmid')[0].files[0];
        var output_rebase = $('#output_rebase')[0].files[0];
        var codon_table = $("#codon_table").val();

        fd.append('or_input_plasmid', or_input_plasmid);
        fd.append('nw_input_plasmid', nw_input_plasmid);
        fd.append('output_rebase', output_rebase);
        fd.append('codon_table', codon_table);

        console.log(codon_table)

        $.ajax({
            type: "POST",
            url: proxy_prefix + "/run_check_output/",
            data: fd,
            processData: false,
            contentType: false,
            dataType: "json",
            success: function (data) {


                if (data['type'] == 'Success') {
                    $('#output_status').remove();

                    $('#loader').remove();

                    $('#input_body').append('Your job id is: '.bold() + data.job_id);

                    $('#input_body').append(components['Check status button']);

                    $('#current_job_id').val(data.job_id);
                    $('#current_job_type').val('check_output');
                } else {
                    alert(data.message);
                    $('#loader').remove();
                }

            },
            error: function (req, status, err) {
                var resp = req.responseJSON;
            }
        });
    } else if (current_page == 'tool_representation') {
        
        var fd = new FormData();

        var input_sequence = $('#input_sequence')[0].files[0];
        var output_rebase = $('#output_rebase')[0].files[0];

        fd.append('input_sequence', input_sequence);
        fd.append('output_rebase', output_rebase);

        $.ajax({
            type: "POST",
            url: proxy_prefix + "/run_tool_representation/",
            data: fd,
            processData: false,
            contentType: false,
            dataType: "json",
            success: function (data) {


                if (data['type'] == 'Success') {
                    $('#output_status').remove();

                    $('#loader').remove();

                    $('#input_body').append('Your job id is: '.bold() + data.job_id);

                    $('#input_body').append(components['Check status button']);

                    $('#current_job_id').val(data.job_id);
                    $('#current_job_type').val('tool_representation');
                } else {
                    alert(data.message);
                    $('#loader').remove();
                }

            },
            error: function (req, status, err) {
                var resp = req.responseJSON;
            }
        });

    } else if (current_page == 'tool_preprocess') {

        var fd = new FormData();

        var input_sequence = $('#input_sequence')[0].files[0];
        var backbone_input_sequence = $('#backbone_input_sequence').val();

        fd.append('input_sequence', input_sequence);
        fd.append('backbone_input_sequence', backbone_input_sequence);

        $.ajax({
            type: "POST",
            url: proxy_prefix + "/run_tool_preprocess/",
            data: fd,
            processData: false,
            contentType: false,
            dataType: "json",
            success: function (data) {

                if (data['type'] == 'Success') {
                    $('#output_status').remove();

                    $('#loader').remove();

                    $('#input_body').append('Your job id is: '.bold() + data.job_id);

                    $('#input_body').append(components['Check status button']);

                    $('#current_job_id').val(data.job_id);
                    $('#current_job_type').val('tool_preprocess');
                } else {
                    alert(data.message);
                    $('#loader').remove();
                }

            },
            error: function (req, status, err) {
                var resp = req.responseJSON;
            }
        });

    } else if (current_page == 'motifsfinder') {

        var bug_ids = JSON.parse($('#bug_ids').html());
        var rebase_ids = JSON.parse($('#rebase_ids').html());

        console.log("Bug ids: ", bug_ids);
        console.log("Rebase ids: ", rebase_ids);

        // Ids:
        // __input_genetic_tool_genus
        // __input_genetic_tool_species
        // __input_genetic_tool_strain
        // __input_genetic_tool_file

        // __input_rebase_genus
        // __input_rebase_species
        // __input_rebase_strain
        // __input_rebase_file

        function callback_Original(bug_ids, rebase_ids) {
            return new Promise((resolve, reject) => {

                var fd = new FormData();

                for (el of bug_ids) {

                    console.log("Bug value:", el);

                    var tmp_gen_tool = $('#' + String(el) + '__input_genetic_tool_file')[0].files[0];
                    fd.append(String(el) + '__input_genetic_tool_file', tmp_gen_tool);

                }

                for (el of rebase_ids) {

                    console.log("REBASE value:", el);

                    var tmp_rebase = $('#' + String(el) + '__input_rebase_file')[0].files[0];
                    fd.append(String(el) + '__input_rebase_file', tmp_rebase);

                }
                if (!!fd.entries().next().value) {
                    console.log('FormData completed!')
                    resolve(fd);
                } else {
                    reject('Empty FormData')
                }
                
            });
        }

        callback_Original(bug_ids, rebase_ids)
            .then(fd => {
                // Loop finished, what to do nexT?
                console.log('Here:', !!fd.entries().next().value)

                // fetch('/run_motifs_finder/', {
                //     method: 'POST',
                //     body: fd
                // }).then(function (response) {
                //     return response.text();
                // }).then(function (text) {
                //     var tmp_text = JSON.parse(text);
                //     console.log(tmp_text);
                // });

                $.ajax({
                    type: "POST",
                    url: proxy_prefix + "/run_motifs_finder/",
                    data: fd,
                    processData: false,
                    contentType: false,
                    dataType: "json",
                    success: function (data) {

                        if (data['type'] == 'Success') {
                            $('#output_status').remove();

                            $('#loader').remove();

                            $('#input_body').append('Your job id is: '.bold() + data.job_id);

                            $('#input_body').append(components['Check status button']);

                            $('#current_job_id').val(data.job_id);
                            $('#current_job_type').val('codon_bias');

                        } else {
                            alert(data.message);
                            $('#loader').remove();
                        }
                    },
                    error: function (req, status, err) {
                        var resp = req.responseJSON;
                    }
                });

            })
            .catch(error => {
                // Error
                console.log(error);
            });

        

        // console.log(!!fd.entries().next().value);

        
    }

}

function check_job_status() {

    var data = new FormData();
    
    var job_id = $("#job_id_check").val();

    if (job_id == undefined) {
        var job_id = $("#current_job_id").val();
    }

    data.append("job_id", job_id);

    fetch(proxy_prefix + '/check_status/', {
        method: 'POST',
        body: data        
    }).then(function (response) {
        return response.text();
    }).then(function (text) {
        var tmp_text = JSON.parse(text);
        if (tmp_text.type == 'Success') {            
            $('#get_job_id_').html(tmp_text['job_status']);
            // alert(tmp_text['job_status'])

            // $('#input_body').append('Your job id is: '.bold() + data.job_id);

            // $('#input_body').append(components['Check status button']);

            $('#check_status_button').remove();

            $('#input_body').append(components['Check status button']);

            $('#curr_output_status').remove();
            // $('#input_body').append('<a id="curr_output_status" style="float: right; padding-right: 10px; padding-top: 10px">' + tmp_text['job_status'] + '</a>');

            // $('#input_body').append('<span class="badge badge-success" id="curr_output_status" style="float: right; padding-right: 30px; padding-top: 20px;">' + tmp_text['job_status'] + '</span>');
            // $('#input_body').append('<span class="badge badge-secondary" id="curr_output_status" style="float: right; margin-right: 10px; width: 100px; height: 50px">' + tmp_text['job_status'] + '</span>');

            if ((tmp_text['job_status'] == 'waiting') | (tmp_text['job_status'] == 'running')) {
                $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #ffbb33; border-color: #FF8800; margin-bottom: 1%; float: right;" id="curr_output_status" disabled>' + tmp_text['job_status'] + '</button>');
            } else if (tmp_text['job_status'] == 'ended successfully') {
                

                $('#check_status_button').remove();

                $('#input_body').append(components['Get output button']);

                $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #00C851; border-color: #007E33; margin-bottom: 1%; float: right;" id="curr_output_status" disabled>' + tmp_text['job_status'] + '</button>');



            } else if (tmp_text['job_status'] == 'failed') {
                $('#input_body').append('<button type="submit" class="btn btn-primary shadow-none" style="background-color: #ff4444; border-color: #CC0000; margin-bottom: 1%; float: right;" id="curr_output_status" disabled>' + tmp_text['job_status'] + '</button>');
            }

            

            

            // <span class="badge badge-success">Success</span>
            // <span class="badge badge-danger">Danger</span>
            // <span class="badge badge-warning">Warning</span>

            


        } else {
            var tmp_text = JSON.parse(text);
            alert(tmp_text.message);
            $('#get_job_id_').html(tmp_text['job_status']);
        }
        
    });
}

function get_output() {
    // $('#main_body').append(components['Loader']);
    // $('#main_body').append(components['Overlay']);
    
    $('#input_body').append(components['Loader']);
    $('#loader').css('float', 'right')

    var fd = new FormData();
    var job_id = $("#job_id_check").val();

    if (job_id == undefined) {
        var job_id = $("#current_job_id").val();
    }

    fd.append('job_id', job_id)

    fetch(proxy_prefix + '/get_output/', {
        method: 'POST',
        body: fd
    }).then(function (resp) {
        return resp.blob();        
    }).then(function (blob) {
        if (blob.type == 'application/zip') {
            download(blob, "out_sytogen_suite.zip"); 
        } else if (blob.type == 'application/json') {
            alert('Incorrect job id');
        }
        
        $('#loader').remove();

    });
}

function change_name(this_obj) {

    var field_name = $(this_obj).attr('id');
    var fileName = $('#' + field_name).val().split("\\");
    var fileName = fileName[fileName.length - 1];
    $('#' + field_name).next('.custom-file-label').html(fileName);

}

function reload() {
    load_current_page();
    $('#remove_backbone').html('false');
    $('#current_considered_slice').empty();
    $('#current_alternative_slice').empty();
    $('#synpl_slice_list').empty();
    $('#configFile').empty();
    $('#selection_object').empty();
    $('#current_job_id').empty();
    $('#current_job_type').empty();
    $('#motifs_dict').empty();
}

function send_message() {
    
    var fd = new FormData();

    var mail = $('#input_email').val();
    var message = $('#comment').val();

    fd.append('mail', mail);
    fd.append('message', message);

    
    fetch(proxy_prefix + '/receive_message/', {
        method: 'POST',
        body: fd
    }).then(function (resp) {
        return resp.json();
    }).then(function (text) {
        
        alert(text["message"]);

        load_current_page();

    });
    
    
    
}

$(document).ready(function () {
    "use strict";

    $("#sidebarToggle").on("click", function (e) {
        e.preventDefault();
        $("body").toggleClass("sb-sidenav-toggled");
    });

    load_current_page();
    
    $("html, body").animate({ scrollTop: 0 }, 10 ^ 10);

});

$(window).on('hashchange', function() {

    load_current_page();

    $("html, body").animate({ scrollTop: 0 }, 10^10);

});

$(window).on('beforeunload', function () {
    $(window).scrollTop(0);
})