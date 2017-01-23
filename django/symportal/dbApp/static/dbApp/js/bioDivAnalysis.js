




    /*  Handle data type selection
    Make the correct upload controls appear */
$("#dataType").change(function(){


    /* if it is fastaQ selected */
    if($(this).val() == 'fastQ'){
        $("#fastQSubmit").show();
        $("#countFileInput").val("")
        $("#fastaFileInput").val("")
        $("#fastQFileInput").val("")
        $("#fastaSubmit").hide();
        $("#fileValidate").hide();

    }

    /* if it is .fasta selected */
    if($(this).val() == ".fasta/count table"){
        $("#fastaSubmit").show();
        $("#countFileInput").val("")
        $("#fastaFileInput").val("")
        $("#fastQFileInput").val("")
        $("#fastQSubmit").hide();
        $("#fileValidate").hide();
    }

    if($(this).val() == ""){
        $("#countFileInput").val("")
        $("#fastaFileInput").val("")
        $("#fastQFileInput").val("")
        $("#fastaSubmit").hide();
        $("#fastQSubmit").hide();
        $("#fileValidate").hide();
    }


});

/* if fastQ file provided show validate button */
$("#fastQFileInput").change(function(){


    if($(this).val() != ""){
        $("#fileValidate").show();
    } else {
        $("#fileValidate").hide();
    }
});

/* if both fasta and count files provided show validate button */
$("#fastaFileInput").change(function(){

    /* do not show the button/hide the button if one of the selction is invalid */
    if($(this).val() != ""){
        if($("#countFileInput").val() != ""){
            $("#fileValidate").show();
        }
    } else {
        $("#fileValidate").hide();
    }

});

$("#countFileInput").change(function(){
    $("#testText").text($(this).val());

    if($(this).val() != ''){
        if($("#fastaFileInput").val() != ''){
        $("#fileValidate").show();
        }
    } else{
        $("#fileValidate").hide();
    }

});