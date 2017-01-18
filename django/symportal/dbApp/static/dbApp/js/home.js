$(document).ready(function () {

    var symPorTxt = "symPortal";
    var symi = 0;

    var symPorTyper = setInterval(function () {
        if (symi != symPorTxt.length) {
            symi += 1;
            $("#symPortalText").text(symPorTxt.substr(0, symi));
        }
        else
        {
         clearInterval(symPorTyper);
         bar()
        }

        console.log(symi);
    }, 50);


    $.when(symPorTyper()).done(function(){



    });

    function bar() {
    $("#hiddenSub").attr("style", "text-align: center; display: block;");
    }

});

/* var subOneTyper = setInterval(function () {
        if (sub1Count != subOneTxt.length) {
            sub1Count += 1;
            $("#sub1").text(text.substr(0, sub1Count));
        }
        else
        {
         clearInterval(subOneTyper);
        }

        console.log(sub1Count);
    }, 50);
    */