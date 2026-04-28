const components = {
    
    "Page title":
        '<h1 class="mt-4" id = "page_title" >__VAR_STR__</h1>',
    "Page subtitle":
        '<ol class="breadcrumb mb-6" id="page_subtitle" style="width:100%;">' +
            '<li class="breadcrumb-item active"></li>' +
            '__VAR_STR__' +
        '</ol>',
    "Image":
        '<img id="page_image" src="'+proxy_prefix+'/static/images/sytogen_suite_pipeline.svg" alt="SyToGenSuite" style="display: block; margin-left: auto; margin-right: auto; width: 90%;">',

    "Documentation":
        '<div class="card">' +
            '<div class="card-header">' +
                'User guide' +
            '</div>' +
            '<div class="card-body">' +
                '<h5 class="card-title">SyToGen Suite</h5>' +
                '<p class="card-text">The SyToGen Suite is a suite of computation tools for genetic engineering.</p>' +
                '<ol id="documentation" style="font-size: 100%;">' +
                    '<li> Build your plasmid(well annotated)</li>' +
                    '<li>Get REBASE output for RM system targets</li>' +
                    '<li>Get your species/strain genome (and in gbk)</li>' +
                    '<li>Run the Syngeneic Tool Estimator</li>' +
                    '<li>Get job id and store it safely</li>' +
                    '<li>Check status of your job when successfully completed</li>' +
                        '<ol type="a">' +
                            '<li>get your output folder</li>' +
                            '<ol type="i">' +
                                '<li>This contains all the generated candidate alternative slces for each slice + an automatically proposed candidate alternative plasmid + codon bias table</li>' +
                            '</ol>' +
                            '<li>produce your syngenic candidate alternvative plasmid through candidate syngenic builder:</li>' +
                            '<ol type="i">' +
                                '<li>Retrieve your output through your job id</li>' +
                                '<li>Choose slices to consider</li>' +
                            '<ol type="1">' +
                                '<li>For the consider slice select the candidate alternative you prefer</li>' +
                                '<li>For all the others set which rank of candidate alternatives to consider</li>' +
                            '</ol>' +
                            '</ol>' +
                        '</ol>' +
                    '<li>Build the primers to produce the syngenic modifications to the sequence</li>' +
                '</ol >' +
            '</div>' + 
        '</div>',

    "Tour card": 
        '<div class="card text-center">' +
            '<div class="card-header">' +
                '<ul class="nav nav-tabs card-header-tabs">' +
                    '<li class="nav-item">' +
                        '<a class="nav-link active" href="#take_a_tour_sytogen_suite" id="take_a_tour_sytogen_suite">SyToGen Suite</a>' +
                    '</li>' +    
                    '<li class="nav-item">' +
                        '<a class="nav-link" href="#take_a_tour_codon_bias_estimator" id="take_a_tour_codon_bias_estimator">Codon bias estimator</a>' +
                    '</li>' +
                    '<li class="nav-item">' +
                        '<a class="nav-link" href="#take_a_tour_synplest" id="take_a_tour_synplest">the Syngeneic Tool Estimator</a>' +
                    '</li>' +
                    '<li class="nav-item">' +
                        '<a class="nav-link" href="#take_a_tour_candidate_syngenic" id="take_a_tour_candidate_syngenic">Candidate syngenic</a>' +
                    '</li>' +
                    '<li class="nav-item">' +
                        '<a class="nav-link" href="#take_a_tour_check_status" id="take_a_tour_check_status">Check status</a>' +
                    '</li>' +
                    '<li class="nav-item">' +
                        '<a class="nav-link" href="#take_a_tour_check_input" id="take_a_tour_check_input">Check input</a>' +
                    '</li>' +
                    '<li class="nav-item">' +
                        '<a class="nav-link" href="#take_a_tour_check_output" id="take_a_tour_check_output">Check output</a>' +
                    '</li>' +
                    '<li class="nav-item">' +
                        '<a class="nav-link" href="#take_a_tour_tool_representation" id="take_a_tour_tool_representation">Tool representation</a>' +
                    '</li>' +
                    '<li class="nav-item">' +
                        '<a class="nav-link" href="#take_a_tour_tool_preprocess" id="take_a_tour_tool_preprocess">Tool preprocess</a>' +
                    '</li>' +
                '</ul>' +
            '</div>' +
            '<div class="card-body">__BODY__</div>' +
        '</div>',
    
    "contacts": 
        '<p>The SyToGen Suite has been developed in collaborations with the <a href="http://segatalab.cibio.unitn.it/people.html" target="_blank">Segata Lab (Unitn)</a> and the <a href="https://johnstonlaboratory.com/contact/" target="_blank">Johnston Lab (Fred Hutch)</a></p>' +

        '<form>' +
            '<div class="form-group">' +
                '<label for="input_email">Email address</label>' +
                '<input type="email" class="form-control" id="input_email" aria-describedby="emailHelp" placeholder="Enter email">' +
                    '<small id="emailHelp" class="form-text text-muted">We\'ll never share your email with anyone else.</small>' +
            '</div>' +
            '<div class="form-group">' +
                '<label for="comment">Message:</label>' +
                '<textarea class="form-control" rows="5" id="comment"></textarea>' +
            '</div>' +
            '<button type="submit" class="btn btn-primary" style="background-color: #337ab7; border-color: #2e6da4;" onclick="send_message()">Send</button>' +
        '</form>',
        
    "Start the tour": 
        '<h1 class="card-title">SyToGen Suite</h1>' +
            '<p class="card-text" style="font-size: 100%;">The SyToGen Suite is a suite of computation tools for genetic engineering.</p>' +
        '<a href="#take_a_tour_codon_bias_estimator" class="btn btn-primary">Start the tour</a>',

    "Tour SyToGen Suite": 
        '<p class="card-text" style="text-align: left; font-size: 100%;">The SyToGen Suite is a suite of computational tool for genetic enginnering. <br></p>',

    'Tour MyMotifs': 
        '<p class="card-text" style="text-align: left; font-size: 100%;">MyMotifs allows to generate a REBASE details file for SyToGen. <br></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Input:</b> <br>' +
            '<ul style="text-align: left; font-size: 100%;">' +
                '<li>REBASE details file</li>' +
                '<li>Other input RM systems</li>' +
            '</ul></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Output:</b> <br>' +
            '<ul style="text-align: left; font-size: 100%;">' +
                '<li>File in REBASE details format</li>' +
            '</ul>' +
        '</p>',
    'Tour MotifsFinder': 
        '<p class="card-text" style="text-align: left; font-size: 100%;">MotifsFinder generates a mapping of RM systems on the input genetic tool. <br></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Input:</b> <br>' +
            '<ul style="text-align: left; font-size: 100%;">' +
                '<li>Multiple genetic tools in GenBank format</li>' +
                '<li>Multiple RM systems files (in REBASE details format)</li>' +
            '</ul></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Output:</b> <br>' +
            '<ul style="text-align: left; font-size: 100%;">' +
                '<li>Map of RM systems on the genetic tools</li>' +
                '<li>Table of number of RM systems in coding and non-coding regions</li>' +
            '</ul>' +
        '</p>',

    "Tour codon bias estimator":
        '<p class="card-text" style="text-align: left; font-size: 100%;">Codon bias estimator returns the codon usage table for the target strain. <br></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Input:</b> <br>' + 
            '<ul style="text-align: left; font-size: 100%;">' +
                '<li>Annotated genome in GenBank format</li>' +
                '<li>Codon table</li>' +
            '</ul></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Output:</b> <br>' +
            '<ul style="text-align: left; font-size: 100%;">' +
                '<li>Codon usage table</li>' +
            '</ul>' +
        '</p>',
    
    "Tour synplest":
        '<p class="card-text" style="text-align: left; font-size: 100%;">The Syngeneic Tool Estimator produces a candidate synegneic alternative for the provided tool.<br></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Input:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
            '<li>Circular genetic tool (plasmid/minicircle) in GenBank format</li>' +
            '<li>REBASE output (for RM systems\' specifications)</li>' +
            '<li>Genome of the target strain (for codon bias estimation) also in GenBank format</li>' +
            '<li>Codon table</li>' +
        '</ul></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Output:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
            '<li>Candidate syngeneic alternative</li>' +
            '<li>Set of .synpl files with ranked evaluated candidate alternative slices, each containing the following information: synonymous change, number of left targets, number of introduced modifications and types of RM systems</li>' +
            '<li>Codon usage table</li>' +
        '</ul>' +
        '</p>',
    
    "Tour candidate syngenic":
        '<p class="card-text" style="text-align: left; font-size: 100%;">Candidate syngenic builder module permits to directly interact with the output candidate alternative slices generated by the Syngeneic Tool Estimator (that can be accessed through the job id) to build a manually-curated syngenic sequence. <br>Steps to follow:</p>' +
        '<ul><li>Retrieve the .synpl file through the job id</li>' + 
        '<li>Click on the slice to consider (from the list on the left)</li>' + 
        '<li>Select a candidate alternative (from the list on the right)</li>' + 
        '<li>When clicking on a different candidate alternative slice the one previously selected is de-selected and will not be considered when building the manually-curated syngenic plasmid</li>' + 
        '<li>For each considered slice no more than one existing candidate alternative can be considered</li>' + 
        '<li>For all the other slices that have not been considered, there is the possibility to consider the original sequence of the slice or the top ranked candidate alternative</li></ul>',
    
    "Tour check status":
        '<p class="card-text" style="text-align: left; font-size: 100%;">Use your job id to check the status of your job (waiting, running, ended successfully or failed). Once the run is successfully completed, get the output folder of your job. You can also do this directly from the check box once it appears after submitting a job.<br></p>' +
        '</p>',

    "Tour check input":
        '<p class="card-text" style="text-align: left; font-size: 100%;">Check input module permits to get information about the input files (in GenBank format) and check whether or not they are correct. <br></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Input:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
        '<li>Circular genetic tool (plasmid/minicircle) or genome (also with multiple sequences) in GenBank format</li>' +
        '</ul></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Output:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
        '<li>Excel file with information about the input sequence and its annotations</li>' +
        '</ul>' +
        '</p>',
    
    "Tour check output":
        '<p class="card-text" style="text-align: left; font-size: 100%;">Check output module is an additional tool for checking that the newly generated syngenic sequence is synonimous to the original one. <br>' +
        'This check is also automatically performed by the Syngeneic Tool Estimator. <br></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Input:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
        '<li>ORIGINAL circular genetic tool (plasmid/minicircle) in GenBank format</li>' +
        '<li>CANDIDATE SYNGENIC circular genetic tool (plasmid/minicircle) in GenBank format</li>' +
        '<li>Codon table code</li>' +
        '</ul></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Output:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
        '<li>Text file with information about comparison between original and candidate syngenic genetic tool</li>' +
        '</ul>' +
        '</p>',

    "Tour tool representation":
        '<p class="card-text" style="text-align: left; font-size: 100%;">Tool representation module permits to have direct access to the tool sequence with all the mapped features and RM target motifs. <br>' +
        'This permits to have a full control on changes of the sequence accounting for all the mapped sequence features.<br></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Input:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
        '<li>Circular genetic tool (plasmid/minicircle) in GenBank format</li>' +
        '<li>REBASE output (for RM target motifs specifications)</li>' +
        '</ul></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Output:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
        '<li>Excel file with the unfolded representation of the plasmid with additional information about mapped target motifs</li>' +
        '</ul>' +
        '</p>',

    "Tour tool preprocess":
        '<p class="card-text" style="text-align: left; font-size: 100%;">Tool preprocess module permits to remove the backbone from the input genetic tool.<br></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Input:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
        '<li>Circular genetic tool (plasmid/minicircle) in GenBank format</li>' +
        '<li>Sequence of the backbone of the genetic tool in Fasta format</li>' +
        '</ul></p>' +
        '<p class="card-text" style="text-align: left; font-size: 100%;"><b>Output:</b> <br>' +
        '<ul style="text-align: left; font-size: 100%;">' +
        '<li>Circular genetic tool without the backbone sequence (minicircle) in GenBank format</li>' +
        '</ul>' +
        '</p>',

    "Input box": 
        '<div class= "card mb-4" id="input_box">' +
            '<div class="card-header">' +
            '<span class= "oi oi-data-transfer-upload" aria-hidden="true"></span>'+
            '</i> Input box<a class="btn btn-light" onclick="reload()"><span class="oi oi-reload"></span></a></div>'  +
            '<div class="card-body" id="input_body">' +
            '</div>' +         
        '</div>',

    "Input file": 
        '<div class="input-group mb-3" id="upload_input_file">' +
            '<div class="input-group-prepend" style="width: 100%">' +
                '<span class="input-group-text">__VAR_LABEL__</span>' +
                '<div class="custom-file">' +
                '<input type="file" class="custom-file-input" id="__VAR_ID__" onchange="change_name(this)">' +
                    '<label class="custom-file-label" for="input_plasmid">Choose file</label>' +
                '</div>' +
            '</div>' +
        '</div>',

    "Add_icon":
        '<div style="margin-bottom: .5%" id="__id_add_icon__"> __info__ <a class="btn btn-light" onclick="__action__"><span class="fas fa-plus"></span></a></div>',

    "Bug_info": 
        '<div class="form-row">' +
            '<div class="col-md-4 mb-3">' +
                '<label for="__bug_genus__">Genus</label>' +
                '<input type="text" class="form-control" id="__bug_genus__" placeholder="Genus" required>' +
                '<div class="invalid-feedback">' +
                    'Please provide a valid genus name.' +
                '</div>' +
            '</div>' +
            '<div class="col-md-4 mb-3">' +
                '<label for="__bug_species__">Species</label>' +
                '<input type="text" class="form-control" id="__bug_species__" placeholder="Species" required>' +
                '<div class="invalid-feedback">' +
                    'Please provide a valid species name.' +
                '</div>' +
            '</div>' +
            '<div class="col-md-4 mb-3">' +
                '<label for="__bug_strain__">Strain</label>' +
                '<input type="text" class="form-control" id="__bug_strain__" placeholder="Strain" required>' +
                '<div class="invalid-feedback">' +
                    'Please provide a valid strain ID.' +
                '</div>' +
            '</div>' +
        '</div>',

    "Remove backbone choice":
        '<form onchange="get_backbone_choice(this)">' + 
            '<div class="form-check-inline">' +
            '    <label class="form-check-label">' +
            '        <input type="radio" class="form-check-input" name="optradio" checked="checked" value="Consider the complete sequence">Consider the complete sequence' +
            '    </label>' +
            '</div>' +
            '<div class="form-check-inline">' +
            '    <label class="form-check-label">' +
            '        <input type="radio" class="form-check-input" name="optradio" value="Remove backbone">Remove backbone' +
            '    </label>' +
        '</form>',

    "MyMotif choice":
        '<form onchange="get_my_motif_choice(this)" id="mymotif_choice">' +
        '<div class="form-check-inline">' +
        '    <label class="form-check-label">' +
        '        <input type="radio" class="form-check-input" name="optradio" checked="checked" value="Modify your REBASE file">Modify your REBASE file' +
        '    </label>' +
        '</div>' +
        '<div class="form-check-inline">' +
        '    <label class="form-check-label">' +
        '        <input type="radio" class="form-check-input" name="optradio" value="Create your REBASE file">Create your REBASE file' +
        '    </label>' +
        '</form>',

    "Caption": '<p id="__CAP_ID__">__VAR_CAPTION__</p>',

    "Backbone sequence": 
        '<div class="input-group" style="padding-bottom: 0.5%" id="__BACKBONE_SEQ__">' + 
            '<div class="input-group-prepend" style="width:100%">' + 
                '<span class="input-group-text">Backbone</span>' + 
                '<textarea class="form-control" aria-label="Backbone sequence"' + 
                    'id="backbone_input_sequence"></textarea>' + 
            '</div>' + 
        '</div>',

    "Input job id":
        '<div class="input-group" style="padding-bottom: 0.5%" id="input_job_id_user">' +
        '<div class="input-group-prepend" style="width:100%">' +
        '<span class="input-group-text">Job id:</span>' +
        '<input type="text" class="form-control" id="__JOB_ID__" placeholder="">' +
        '</div>' +
        '</div>',
    
    "Submit caption": '<p>Note that it may take some time for files to load. Wait a few minutes after running your job.</p>',
    "Submit button": '<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4;" id="submit_button" onclick="submit_job()">Submit</button>',
    "Retrieve folder button": "",
    "Check status button": "",
    "Get output button": "",
    "Generated slices": "",
    "Candidate alternatives": "",
    "Superuser": "",

    "Loader": '<button class="btn btn-primary" type="button" id="loader" style="background-color: #337ab7; border-color: #2e6da4" disabled><span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span>Loading...</button>',

    "Check status button": '<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4; margin-bottom: 1%; float: right;" id="check_status_button" onclick="check_job_status()">Check status</button>',

    "Get output button": '<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4; margin-left: 0.1%; margin-bottom: 1%; float: right;" id="get_output_button" onclick="get_output()">Get output</button>',

    "Retrieve generated slices button": '<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4; margin-left: 0.1%; margin-bottom: 1%; margin-right: 0.5%" id="retrieve_list_of_slices_button" onclick="retrieve_list_of_slices()">Submit</button>',

    "Codon table choice": 
        '<div class="input-group mb-3">' +
            '<div class="input-group-prepend" style="width: 100%">' +
                '<label class="input-group-text" for="codon_table">Codon table</label>' +
                '<select class="custom-select" id="codon_table">' +
                    '<option value="1" selected>1. Standard Code</option>' +
                    '<option value="2">2. Vertebrate Mitochondrial Code</option>' +
                    '<option value="3">3. Yeast Mitochondrial Code</option>' +
                    '<option value="4">4. Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code</option>' +
                    '<option value="5">5. Invertebrate Mitochondrial Code</option>' +
                    '<option value="6">6. Ciliate, Dasycladacean and Hexamita Nuclear Code</option>' +
                    '<option value="9">9. Echinoderm and Flatworm Mitochondrial Code</option>' +
                    '<option value="10">10. Euplotid Nuclear Code</option>' +
                    '<option value="11">11. Bacterial, Archaeal and Plant Plastid Code</option>' +
                    '<option value="12">12. Alternative Yeast Nuclear Code</option>' +
                    '<option value="13">13. Ascidian Mitochondrial Code</option>' +
                    '<option value="14">14. Alternative Flatworm Mitochondrial Code</option>' +
                    '<option value="16">16. Chlorophycean Mitochondrial Code</option>' +
                    '<option value="21">21. Trematode Mitochondrial Code</option>' +
                    '<option value="22">22. Scenedesmus obliquus Mitochondrial Code</option>' +
                    '<option value="23">23. Thraustochytrium Mitochondrial Code</option>' +
                    '<option value="24">24. Rhabdopleuridae Mitochondrial Code</option>' +
                    '<option value="25">25. Candidate Division SR1 and Gracilibacteria Code</option>' +
                    '<option value="26">26. Pachysolen tannophilus Nuclear Code</option>' +
                    '<option value="27">27. Karyorelict Nuclear Code</option>' +
                    '<option value="28">28. Condylostoma Nuclear Code</option>' +
                    '<option value="29">29. Mesodinium Nuclear Code</option>' +
                    '<option value="30">30. Peritrich Nuclear Code</option>' +
                    '<option value="31">31. Blastocrithidia Nuclear Code</option>' +
                    '<option value="33">33. Cephalodiscidae Mitochondrial UAA-Tyr Code</option>' +
                '</select>' +
            '</div>' +
        '</div>',

    'Output status':
        '<div id="output_status">' +
            '<div class="card mb-4">' +
                '<div class="card-header"><i class="fas fa-chart-bar mr-1"></i>Job status</div>' +
                '<div class="alert alert-success" role="alert">' +
                '<h4 class="alert-heading">Your job is:</h4>' +
                '<div id="get_job_id_"></div>' +
                '</div>' +
                '<div class="form-row" style="width: 99%; margin-left: 0.5%">' +
                '<div class="form-group col-xl-12">' +
                '<input type="text" class="form-control" id="__JOB_ID__" placeholder="">' +
                '</div>' +
                '</div>' +
                '<div class="btn-group" role="group" style="width: 99%; margin-left: 0.5%">' +
                '<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4; margin-bottom: 1%;" id="__BUTTON_ID__" onclick="check_job_status()">Check status</button>' +
                '<button type="submit" class="btn btn-primary shadow-none" style="background-color: #337ab7; border-color: #2e6da4; margin-left: 0.1%; margin-bottom: 1%;" id="__BUTTON_ID__" onclick="get_output()">Get output</button>' +
            '</div>' +
        '</div>',

    'Candidate syngenic input': '<div class="input-group-prepend">' +
        '<span class="input-group-text" style="margin-bottom: 0.5%">Job id: </span>' +
        '<input type="text" class="form-control" id="syn_cand_builder_job_id" placeholder="" aria-label="Username" aria-describedby="basic-addon1" style="margin-bottom: 0.5%">',
    
    'Candidate syngenic':

        '<div class="row" style="margin-top: 1%">' +

            '<div id="list_of_synpl_files_box" class="col-xl-6" style="visibility: hidden">' +
                '<div class="card mb-4">' +
                    '<div class="card-header"><i class="fas fa-table"></i>List of generated slices (.synpl files)</div>' +
                        '<div class="card-body">' +
                            '<div class="panel-body">' +
                            '<div class="list-group" id="select_slice" style="font-family:monospace; max-height: 450px; overflow-y:scroll;"></div>' +
                        '</div>' +
                    '</div>' +
                '</div>' +
            '</div>' +
            
            '<div id="list_of_alternative_files_box" class="col-xl-6" style="visibility: hidden">' +
                '<div class="card mb-4">' +
                    '<div class="card-header"><i class="fas fa-table"></i>List of the candidate syngenic alternative for the selected slice</div>' +
                    '<div class="card-body">' +
                    '<div class="panel-body">' +
                        '<div class="list-group" id="temporary_select_slice" style="font-family:monospace; height: 450px; overflow-y:scroll;"></div>' +
                    '</div>' +
                    '</div>' +
                '</div>' +
            '</div>' +
            
        '</div>',

    'Candidate syngenic oth input': 
        '<div class="input-group mb-3" style="margin: 0px 0px 200px 0px" id="other_alternatives">' +
            '<div class="input-group-prepend">' +
                '<span class="input-group-text">All the rest candidate alternatives</span>' +
            '</div>' +
            '<select class="custom-select" id="other_alternatives_input" placeholder="Positive strand (forward: +) 5\' - 3\'">' +
                '<option value="0">Original sequence</option>' +
                '<option value="1">Top ranked (best) candidate alternative</option>' +
            '</select>' +
        '</div>',
    'Candidate syngenic oth button': '<button id="build_candidate_syngenic_button_" style="background-color: #337ab7; border-color: #2e6da4; margin-left: 0.1%; margin-bottom: 1%; margin-top: 0.5%; margin-right: 0.5%" type="submit" class="btn btn-primary shadow-none" method="POST" onclick="build_candidate_syngenic()">Build the candidate SynGenic</button>'
        
}

camth = "rec_seq	Methyl-Modification introduced (5'-3')	meth_base	comp_meth_base	comp_meth_type	meth_type	availability\n" +
"ATCGAT	ATCGm6AT	5	2	6	6	M.ClaI(Takara Bio); M.BseCI(Minotech Biotechnology)\n" +
"AAGCTT	m6AAGCTT	1	6	6	6	M.HindIII(Takara Bio)\n" +
"GGATCC	GGATm4CC	5	2	4	4	M.BamHI(Takara Bio); M.BamHI(New England Biolabs)\n" +
"GAATTC	GAm6ATTC	3	4	6	6	M.EcoRI(Takara Bio); M.EcoRI(New England Biolabs); M.EcoRI(Nippon Gene);\n" +
"GGATG	GGm6ATG	3	3	6	6	M3.BstF5I(SibEnzyme)\n" +
"GCNGC	Gm5CNGC	2	4	5	5	M.Fsp4HI(SibEnzyme)\n" +
"AGCT	AGm5CT	3	2	5	5	M.AluI(Takara Bio); M.AluI(New England Biolabs)\n" +
"TCGA	TCGm6A	4	1	6	6	M.TaqI(New England Biolabs)\n" +
"GATC	Gm6ATC	2	3	6	6	M.EcoKDam(New England Biolabs)\n" +
"GGCC	GGm5CC	3	2	5	5	M.HaeIII(Takara Bio); M.HaeIII(New England Biolabs)\n" +
"GCGC	Gm5CGC	2	3	5	5	M.HspAI(SibEnzyme)\n" +
"CCGG	m5CCGG	1	4	5	5	M.MspI(New England Biolabs)\n" +
"CCGG	Cm5CGG	2	3	5	5	M.HpaII(Takara Bio); M.HpaII(New England Biolabs)\n" +
"GC	Gm5C	2	1	5	5	M.CviPI(New England Biolabs)\n" +
"CG	m5CG	1	2	5	5	M.SssI(New England Biolabs); M.SssI(ThermoFischer); M.SssI(Zymo Research)\n" +
"A	m6A	1	1	6	6	M.EcoGII(New England Biolabs)"