$(document).on('click', '.gene-name', function() {
  var geneName = $(this).text();


  // Requête AJAX pour obtenir l'ID Ensembl
  $.getJSON('https://www.proteinatlas.org/api/search_download.php', {
        search: geneName,
        format: "json",
        columns: "g,eg",
        compress: "no"
  }, function(data) {
        if(data && data.length > 0 && data[0].Ensembl) {
          var ensemblId = data[0].Ensembl;
          var url = "https://www.proteinatlas.org/" + ensemblId;
          window.open(url, '_blank');
        } else {
          alert("Gene not found on Protein Atlas.");
        }
  });
});