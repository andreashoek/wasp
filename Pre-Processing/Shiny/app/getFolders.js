     document.getElementById("fileIn").addEventListener("change", function(e) {

            let files = e.target.files;
            var arr = new Array(files.length*2);
            for (let i=0; i<files.length; i++) {

			arr[i] = files[i].webkitRelativePath;
            arr[i+files.length] = files[i].name;


            }

            Shiny.onInputChange("mydata", arr);

    });