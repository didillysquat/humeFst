$("#testText").text("Happy");


$('#seqModal').on('show.bs.modal', function (event) {
var button = $(event.relatedTarget); // Button that triggered the modal
var seqId = button.data('seq_id'); // Extract info from data-* attributes
var modal = $(this);
//modal.find('.modal-title').text(seqId);
$("#testText").html(seqId);
    // Perople suggested putting the argument for dataType in the get method but that was making it difficult
    // to access what was being passes from python
    // I have now left this out and so when I pass the object it comes as a string
    // I will use JSON.parse(data) to convert this into an accessible onject as seen in the large comment at the bottom
    $.get('seqModal/', {seq_id:seqId},  function(data){

            var json = JSON.parse(data);
            modal.find('.modal-title').text(json.title);
            modal.find('.fasta-title').text(json.fastaTitle);
            modal.find('.fasta-sequence').text(json.sequence);
      });


});









/*
Below is a JSON string.

JSON String
{
    "name": "mkyong",
    "age": 30,
    "address": {
        "streetAddress": "88 8nd Street",
        "city": "New York"
    },
    "phoneNumber": [
        {
            "type": "home",
            "number": "111 111-1111"
        },
        {
            "type": "fax",
            "number": "222 222-2222"
        }
    ]
}
To access the JSON object in JavaScript, parse it with JSON.parse(), and access it via “.” or “[]”.

JavaScript
<script>
       var data = '{"name": "mkyong","age": 30,"address": {"streetAddress": "88 8nd Street","city": "New York"},"phoneNumber": [{"type": "home","number": "111 111-1111"},{"type": "fax","number": "222 222-2222"}]}';

	var json = JSON.parse(data);

	alert(json["name"]); //mkyong
	alert(json.name); //mkyong

	alert(json.address.streetAddress); //88 8nd Street
	alert(json["address"].city); //New York

	alert(json.phoneNumber[0].number); //111 111-1111
	alert(json.phoneNumber[1].type); //fax

	alert(json.phoneNumber.number); //undefined
</script>
*/
















