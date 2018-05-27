

//Some global variables
var graph = 0; // to know if a feature graph has been made in the actual session (0 = not, 1 = yes)

//EVENT LISTENERS
//Input tab event listener
document.getElementById('input_tab').addEventListener("click", function(){
   document.getElementById('File_input_section').style.display = "block";
   document.getElementById('FM_section').style.display = "none";
   document.getElementById('Export_section').style.display = "none";
});

//Feature map tab event listener
document.getElementById('FM_tab').addEventListener("click", function(){
   document.getElementById('File_input_section').style.display = "none";
   document.getElementById('Export_section').style.display = "none";
   document.getElementById('FM_section').style.display = "block";
});

//Export tab event listener
document.getElementById('Export_tab').addEventListener("click", function(){
   document.getElementById('File_input_section').style.display = "none";
   document.getElementById('FM_section').style.display = "none";
   document.getElementById('Export_section').style.display = "block";
});

//Advanced options button event listener
var advanced_options_button = 0;
document.getElementById('Advanced_options_button').addEventListener("click", function(){
   if(advanced_options_button == 0){
      document.getElementById('Advanced_options_div').style.display = "block";
      document.getElementById('Advanced_options_button').innerHTML = "Hide advanced options";
      advanced_options_button = 1;
   }
   else{
      document.getElementById('Advanced_options_div').style.display = "none";
      document.getElementById('Advanced_options_button').innerHTML = "View advanced options";
      advanced_options_button = 0;
   }
});

// Go button event listener
document.getElementById("Go_button").addEventListener("click",function(){
   document.getElementById("File_input_section").style.display = "none";
   draw_feature_map();
   document.getElementById("FM_section").style.display = "block";

   graph = 1;

   document.getElementById("input_tab").classList.remove("current");
   document.getElementById("FM_tab").classList.add("current");
});

// Reset button event listener
document.getElementById("Reset_button").addEventListener("click",function(){
   reset_fields();
});

//CODE ON PROGRESS --- START  --- File uploading option
// var f_name;
// function handleFileSelect(evt) {
// var files = evt.target.files; // FileList object
//
// // files is a FileList of File objects. List some properties.
//
// var output = [];
// for (var i = 0, f; f = files[i]; i++) {
//    output.push('<li><strong>', escape(f.name), '</strong> (', f.type || 'n/a', ') - ',
//                f.size, ' bytes, last modified: ',
//                f.lastModifiedDate ? f.lastModifiedDate.toLocaleDateString() : 'n/a',
//                '</li>');
//
//    console.log("NOMBRE DE ARCHIVO: "+ f.name);
// }
// document.getElementById('list').innerHTML = '<ul>' + output.join('') + '</ul>';
//
// }
//
// document.getElementById('files').addEventListener('change', handleFileSelect, false);

//CODE ON PROGRESS --- END --- File uploading option

//JAVASCRIPT FUNCTIONS
//Function to calculate the with of a text string according to a given font and size
function calculate_text_width(text, font, size){
   var c = document.getElementById("myCanvas");
   var ctx = c.getContext("2d");
   ctx.font = size+"px "+font;
   var txt = text;
   var size = ctx.measureText(txt).width;
   return size;
}

//Function to be called by the reset button to restore all the values of the input tab (basic and advanced options)
function reset_fields(){
   //Reset of basic options
   //Reset title of graph field
   document.getElementById("user_title").value = "";

   //Reset Upload your file field
   document.getElementById("FileUpload_text").value = "";

   //Reset Select a format field
   document.getElementById("Select_format").value = "feature map";

   //Reset Feature list field
   document.getElementById("Feature_list").value = "";

   //Reset of advanced options
   //Reset of title position
   document.getElementById("title_left").checked = "true";

   //Reset of legends view labels
   document.getElementById("view_legend").checked = "true";

   //Reset of legends position
   document.getElementById("legend_top").checked = "true";

   //Reset of view Scalebar
   document.getElementById("view_scale").checked = "true";

   //Reset of view score
   document.getElementById("view_score").checked = "true";

}

// Function to draw the feature map svg
function draw_feature_map(){

   //Get values of all the basic and advanced options to make the first iteration of the feature map graph
   var user_title = document.getElementById("user_title").value;
   var json_file = document.getElementById('FileUpload').files[0]; // File path

   var title_position_left = document.getElementById('title_left');
   var title_position_center = document.getElementById('title_center');
   var title_position_right = document.getElementById('title_right');

   //Advanced options variables
   var view_legend = document.getElementById('view_legend');
   var legend_top = document.getElementById('legend_top');
   var legend_bottom = document.getElementById('legend_bottom');
   var legend_left = document.getElementById('legend_left');
   var legend_right = document.getElementById('legend_right');

   var view_scale = document.getElementById('view_scale')
   var view_score = document.getElementById('view_score');


   //Opening of the file, important to do it first because some calculations need to be performed before the creation of the SVG
   d3.json("js/prueba_3.json", function(data){

      //Printing coordinates for the SVG graph
      var p_title_x = 25;
      var p_title_y = 20;
      var p_title_width = 0;
      var p_title_height = 0;
      var p_title_margin_left = 0;
      var p_title_space = 20;

      var p_legends_top_x = 0;
      var p_legends_top_y = 0;
      var p_legends_top_width = 0;
      var p_legends_top_height = 0;

      var p_legends_right_x = 0;
      var p_legends_right_y = 0;
      var p_legends_right_width = 0;
      var p_legends_right_y_height = 0;

      var p_legends_left_x = p_title_x;
      var p_legends_left_y = 0;
      var p_legends_left_width = 0;
      var p_legends_left_height = 0;

      var p_legends_bottom_x = 0;
      var p_legends_bottom_y = 0;
      var p_legends_bottom_width = 0;
      var p_legends_bottom_height = 0;

      var p_graph_x = 0;
      var p_graph_y = 0;
      var p_graph_width = 0;
      var p_graph_height = 0;

      var total_svg_width = document.getElementById('feature_maps').offsetWidth;  // According to the size of the div, this depends of the grid configuration on the HTML file

      //Calculate the height  of the SVG according to number of genes in the json file
      var gene_num = 0;
      var dedicated_space = 100;
      var title_space = 100;

      for(i=0; i<data.length; i++){
         gene_num = gene_num +1;
      }

      //Calculate height according to file and user preferences
      var total_svg_height = (gene_num * dedicated_space + title_space);

      //Creation of the svg
      var graph_svg = d3.select("#feature_maps")
         .append("svg")
         .attr("x",0)
         .attr("y",0)
         .attr("width",total_svg_width)
         .attr("height",total_svg_height);

      // TITLE BAR CODE
      var g_title_x;
      var g_title_y;
      var g_title_anchor;

      //If title position LEFT
      if (title_left.checked == true){
         g_title_x = p_title_x + p_title_margin_left;
         g_title_y = p_title_y;
         g_title_anchor = "start";
      }


      //If title position is CENTER
      if (title_center.checked == true){
         g_title_x = total_svg_width/2;
         g_title_y = p_title_y;
         g_title_anchor = "middle";
      }

      //If title position id RIGHT
      if (title_right.checked == true){
         g_title_x = total_svg_width -20;
         g_title_y = p_title_y;
         g_title_anchor = "end";
      }

      //Print of the application title
      graph_svg.append("text")
         .attr("x",g_title_x)
         .attr("y",g_title_y)
         .attr("font-size", 16)
         .attr("text-anchor",g_title_anchor)
         .attr("font-family", "Lato")
         .attr("fill","rgb(51,51,51)")
         .text("RSAT Feature map");

      //User title positioning
      var title_user_x = g_title_x;
      var title_user_y = g_title_y + p_title_space;



      //Print of the User title
      graph_svg.append("text")
         .attr("x",title_user_x)
         .attr("y",title_user_y)
         .attr("font-size", 16)
         .attr("text-anchor", g_title_anchor)
         .attr("font-family", "Lato")
         .attr("fill","rgb(51,51,51)")
         .text(user_title);

      //LEGENDS BAR
      //Get all the legends and extract discard the repeated ones
      var label_array_total = [];

      for(i=0; i < data.length; i++){
         for(f=0; f< data[i].Features.length; f++){
            label_array_total.push(data[i].Features[f].seq);
            //console.log(data[i].Features[f].seq);
         }
      }
      var sorted_label_array_total = label_array_total.sort();

      //Color asignment RANDOM
      var label_array_color = [];// Array of labels and randonmly asigned colors.
      var label_color;
      var rand_color;

      for(i=0; i < sorted_label_array_total.length; i++){
         if(sorted_label_array_total[i] != sorted_label_array_total[i+1]){
            rand_color = randomColor();
            label_color = {label:sorted_label_array_total[i], color:rand_color};
            label_array_color.push(label_color);
            //console.log("L = "+ label_color.label + " C = "+ label_color.color);
         }
      }

      //Legends bar positioning

      //If view legends is checked
      if(view_legend.checked == true){
         //IF LABELS ON TOP
         if(legend_top.checked == true){
            var legends_x = 25;
            var legends_y = 80;

            //Legends title
               graph_svg.append("text")
                  .attr("x",legends_x)
                  .attr("y",legends_y)
                  .attr("font-size", 16)
                  .attr("text-anchor","start")
                  .attr("font-family", "Lato")
                  .attr("fill","rgb(51,51,51)")
                  .text("Legends");

               //Code for printing the legends automatically

               var max_size = 0;
               var rect_size = 10;

               //Get the size of the longest label
               for(i=0; i< label_array_color.length; i++){
                  size_2 = calculate_text_width(label_array_color[i].label,"Lato",12);

                  if (max_size < size_2){
                     max_size = size_2;
                  }
               }
               //console.log("Size " + max_size);

               var tab_size = 40;
               var auto_columns = Math.floor((total_svg_width / (tab_size + max_size + rect_size)));

               //console.log("Auto_columns "+ auto_columns);

               //Printing of the labels and their color squares
               var labels_x = legends_x;
               var labels_y = legends_y + 20;

               var labels_x_jump = rect_size + tab_size + max_size;
               var labels_y_jump = 20;

               var cont = 0;

               for(i=0; i< label_array_color.length; i++){
                  graph_svg.append("rect")
                     .attr("x",labels_x)
                     .attr("y",labels_y)
                     .attr("width",rect_size)
                     .attr("height",rect_size)
                     .attr("fill",label_array_color[i].color);

                  graph_svg.append("text")
                     .attr("x",labels_x + rect_size +4)
                     .attr("y",labels_y + rect_size -1)
                     .attr("font-size", 12)
                     .attr("text-anchor","start")
                     .attr("font-family", "Lato")
                     .attr("fill","rgb(51,51,51)")
                     .text(label_array_color[i].label);

                  cont = cont +1;
                  labels_x = labels_x + labels_x_jump;

                  if(cont >= auto_columns){
                     labels_y = labels_y + 20;
                     labels_x = legends_x;
                     cont = 0;
                  }
               }

         }
         // IF LABELS ON THE right
         if(legend_right.checked == true){

            var max_size = 0;
            var rect_size = 10;

            //Get the size of the longest label
            for(i=0; i< label_array_color.length; i++){
               size_2 = calculate_text_width(label_array_color[i].label,"Lato",12);

               if (max_size < size_2){
                  max_size = size_2;
               }
            }

            var legends_x = (total_svg_width - 20 - max_size); // NEED TO BE CORRECTED
            var legends_y = 80;

            //Legends title
            graph_svg.append("text")
               .attr("x",legends_x)
               .attr("y",legends_y)
               .attr("font-size", 16)
               .attr("text-anchor","start")
               .attr("font-family", "Lato")
               .attr("fill","rgb(51,51,51)")
               .text("Legends");

               var labels_y = legends_y + 20;
               var labels_y_jump = 20;
               var labels_x = legends_x;

               for(i=0; i< label_array_color.length; i++){
                  graph_svg.append("rect")
                     .attr("x",labels_x)
                     .attr("y",labels_y)
                     .attr("width",rect_size)
                     .attr("height",rect_size)
                     .attr("fill",label_array_color[i].color);

                  graph_svg.append("text")
                     .attr("x",labels_x + rect_size +4)
                     .attr("y",labels_y + rect_size -1)
                     .attr("font-size", 12)
                     .attr("text-anchor","start")
                     .attr("font-family", "Lato")
                     .attr("fill","rgb(51,51,51)")
                     .text(label_array_color[i].label);


                  labels_y = labels_y + labels_y_jump;

               }

         }
         // IF LABELS ON THE left
         if(legend_left.checked == true){

            var max_size = 0;
            var rect_size = 10;

            //Get the size of the longest label
            for(i=0; i< label_array_color.length; i++){
               size_2 = calculate_text_width(label_array_color[i].label,"Lato",12);

               if (max_size < size_2){
                  max_size = size_2;
               }
            }

            var legends_x = p_legends_left_x; // NEED TO BE CORRECTED
            var legends_y = 80;

            //Legends title
            graph_svg.append("text")
               .attr("x",legends_x)
               .attr("y",legends_y)
               .attr("font-size", 16)
               .attr("text-anchor","start")
               .attr("font-family", "Lato")
               .attr("fill","rgb(51,51,51)")
               .text("Legends");

               var labels_y = legends_y + 20;
               var labels_y_jump = 20;
               var labels_x = legends_x;

               for(i=0; i< label_array_color.length; i++){
                  graph_svg.append("rect")
                     .attr("x",labels_x)
                     .attr("y",labels_y)
                     .attr("width",rect_size)
                     .attr("height",rect_size)
                     .attr("fill",label_array_color[i].color);

                  graph_svg.append("text")
                     .attr("x",labels_x + rect_size +4)
                     .attr("y",labels_y + rect_size -1)
                     .attr("font-size", 12)
                     .attr("text-anchor","start")
                     .attr("font-family", "Lato")
                     .attr("fill","rgb(51,51,51)")
                     .text(label_array_color[i].label);


                  labels_y = labels_y + labels_y_jump;

               }

         }

      // IF LABELS AT THE BOTTOM

      }


   //    //Legends positioning
   //    var legends_x = 25;
   //    var legends_y = 80;
   //
   //    //Legends title
   //    graph_svg.append("text")
   //       .attr("x",legends_x)
   //       .attr("y",legends_y)
   //       .attr("font-size", 16)
   //       .attr("text-anchor","start")
   //       .attr("font-family", "Lato")
   //       .attr("fill","rgb(51,51,51)")
   //       .text("Legends");
   //
   //    //Code for printing the legends automatically
   //
   //    var max_size = 0;
   //    var rect_size = 10;
   //
   //    //Get the size of the longest label
   //    for(i=0; i< label_array_color.length; i++){
   //       size_2 = calculate_text_width(label_array_color[i].label,"Lato",12);
   //
   //       if (max_size < size_2){
   //          max_size = size_2;
   //       }
   //    }
   //    console.log("Size " + max_size);
   //
   //    var tab_size = 40;
   //    var auto_columns = Math.floor((svg_width / (tab_size + max_size + rect_size)));
   //
   //    console.log("Auto_columns "+ auto_columns);
   //
   //    //Printing of the labels and their color squares
   //    var labels_x = legends_x ;
   //    var labels_y = legends_y + 20;
   //
   //    var labels_x_jump = rect_size + tab_size + max_size;
   //    var labels_y_jump = 20;
   //
   //    var cont = 0;
   //
   //    for(i=0; i< label_array_color.length; i++){
   //       graph_svg.append("rect")
   //          .attr("x",labels_x)
   //          .attr("y",labels_y)
   //          .attr("width",rect_size)
   //          .attr("height",rect_size)
   //          .attr("fill",label_array_color[i].color);
   //
   //       graph_svg.append("text")
   //          .attr("x",labels_x + rect_size +4)
   //          .attr("y",labels_y + rect_size -1)
   //          .attr("font-size", 12)
   //          .attr("text-anchor","start")
   //          .attr("font-family", "Lato")
   //          .attr("fill","rgb(51,51,51)")
   //          .text(label_array_color[i].label);
   //
   //       cont = cont +1;
   //       labels_x = labels_x + labels_x_jump;
   //
   //       if(cont >= auto_columns){
   //          labels_y = labels_y + 20;
   //          labels_x = legends_x;
   //          cont = 0;
   //       }
   //    }
   //
   //
   //Print of the main scale bar
   var scale_bar_x = 120;
   var scale_bar_y = labels_y + 50;
   var scale_width = total_svg_width - 20 - scale_bar_x;
   //
   graph_svg.append("line")
      .attr("x1",scale_bar_x)
      .attr("y1",scale_bar_y)
      .attr("x2",total_svg_width -50)
      .attr("y2",scale_bar_y)
      .attr("stroke","#000000")
      .attr("stroke-width",.5);

   // //Obtener el mayor rango de los genes para calcular el rango de la main scale_bar
    var range = 0;
   //
    for(i=0; i<data.length; i++){
       if(data[i].Start > range){
          range = data[i].Start;
       }
    }
   // console.log("range "+ range);
   //
   // //Graficar ticks de la barra principal
   // //Primero hacer calculo de la escala
    var line_width = (total_svg_width - scale_bar_x - 50);
    var tick_space = 50; // Puede ser definida por el usuario
    var tick_jump = ((tick_space * line_width) / range);
    var value = range;
   //
    while(scale_bar_x < (total_svg_width - 49)){
   //
       graph_svg.append("line")
          .attr("x1",scale_bar_x)
          .attr("y1",scale_bar_y-3)
          .attr("x2",scale_bar_x)
          .attr("y2",scale_bar_y+3)
          .attr("stroke","#000000")
          .attr("stroke-width",1);
   //
       graph_svg.append("text")
         .attr("x",scale_bar_x)
         .attr("y",scale_bar_y+12)
         .attr("font-size", 10)
         .attr("text-anchor","middle")
         .attr("font-family", "Lato")
         .attr("fill","rgb(51,51,51)")
         .text(value);
   //
       scale_bar_x = scale_bar_x + tick_jump;
       value = value - tick_space;
    }
   //
   // //Print of each of the gene names
   var gene_name_x;
   var gene_name_y = p_title_y + 150;
   var gene_jump = 80;
   //
   graph_svg.selectAll("text.Gene")
      .data(data)
      .enter()
      .append("text")
      .attr("x",p_title_x)
      .attr("y", function(d,i){
         return gene_name_y + ((i+1)*gene_jump);
      })
      .attr("font-size", 16)
      .attr("text-anchor","start")
      .attr("font-family", "Lato")
      .attr("fill","rgb(51,51,51)")
      .text(function(d){
         return d.Gene;
      });
   //
   // //Print each of the gene horizontal lines
   var gene_line_x2 = total_svg_width - 50;
   var gene_line_y ;
   var start_y = 95;
   //
   graph_svg.selectAll("line.bar")
      .data(data)
      .enter()
      .append("line")
      .attr("x1", total_svg_width - 50)
      .attr("y1",function(d,i){
         return (gene_name_y + ((i+1)*gene_jump))-10;
         })
      .attr("x2",function(d){
         return ((total_svg_width - 50) - ((d.Start * line_width)/range));
         })
      .attr("y2",function(d,i){
         return (gene_name_y + ((i+1)*gene_jump))-10;
         })
      .attr("stroke","#333333")
      .attr("stroke-width",1);

   // //Print ticks on each line
   //
   var gene_line_ticks_y;
   var gene_line_ticks_x = total_svg_width - 50;

   for(i=0; i < data.length; i ++){

      gene_line_ticks_y = gene_name_y + ((i+1)*gene_jump)-10;

      //Tick on cero
       graph_svg.append("line")
          .attr("x1",total_svg_width - 50)
          .attr("y1",gene_line_ticks_y - 3)
          .attr("x2",total_svg_width - 50)
          .attr("y2",gene_line_ticks_y + 3)
          .attr("stroke","#333333")
          .attr("stroke-width",1);

      //Tick on Start
      graph_svg.append("line")
        .attr("x1",total_svg_width - 50 - ((data[i].Start * line_width)/range))
        .attr("y1",gene_line_ticks_y - 3)
        .attr("x2",total_svg_width - 50 - ((data[i].Start * line_width)/range))
        .attr("y2",gene_line_ticks_y + 3)
        .attr("stroke","#333333")
        .attr("stroke-width",.5);

      while(gene_line_ticks_x > (total_svg_width - 50 - ((data[i].Start * line_width)/range))){
         graph_svg.append("line")
             .attr("x1",gene_line_ticks_x)
             .attr("y1",gene_line_ticks_y - 3)
             .attr("x2",gene_line_ticks_x)
             .attr("y2",gene_line_ticks_y + 3)
             .attr("stroke","#333333")
             .attr("stroke-width",.5);

             gene_line_ticks_x = gene_line_ticks_x - ((50 * line_width)/range);
      }
      gene_line_ticks_x = total_svg_width - 50;
   }

//Print each of the features

   for(i=0; i < data.length; i++){
      for(f=0; f < data[i].Features.length; f++){
         graph_svg.append("rect")
            .attr("x",((total_svg_width - 50) - (data[i].Features[f].start * (line_width/range))))
            .attr("y",function(){
               if(data[i].Features[f].strand == "DR"){
                  return (gene_name_y + ((i+1)*gene_jump) - (data[i].Features[f].score * 3)-10);
               }
               else{
                  return (gene_name_y + ((i+1)*gene_jump)-10);
               }
            })
            .attr("width",(((data[i].Features[f].start - data[i].Features[f].end) * (line_width/range))))
            .attr("height", (data[i].Features[f].score * 3))
            .attr("fill", function(){
               for(c=0; c < label_array_color.length; c++){

                  if(data[i].Features[f].seq == label_array_color[c].label){
                     console.log("SEC: " + data[i].Features[f].seq + "LABEL: " + label_array_color[c].label + "Color " + label_array_color[c].color);
                     return label_array_color[c].color;
                  }
               }
            });
      }
   }
      //Code for graphing features using D3.JS -- Not working
      // graph_svg.selectAll("rect.features")
      //    .data(data[i].Features)
      //    .enter()
      //    .append("rect")
      //    .attr("x",function(d){
      //       return ((total_svg_width - 50) - ((d.start* line_width)/range));
      //    })
      //    .attr("y",function(d,i){
      //       return (gene_name_y + ((i+1)*gene_jump))-10;
      //       })
      //    .attr("width",function(d){
      //       return (10);
      //    })
      //    .attr("height",rect_size)
      //    .attr("fill","blue");


   // graph_svg.selectAll("rect.features")
   //    .data(data)
   //    .enter()
   //    .append("rect")
   //    .attr("x",function(d){
   //       return d.Features.start + 700;
   //    })
   //    .attr("y",function(d,i){
   //       return (gene_name_y + ((i+1)*gene_jump))-10;
   //       })
   //    .attr("width",rect_size)
   //    .attr("height",rect_size)
   //    .attr("fill","blue");

});//closing of the opening of the file option
}

//JQuery functions
//Function for the top menu (input, feature map, export)
$('.breadcrumb-counter-nav-item').click(function () {
  $('.breadcrumb-counter-nav-item').removeClass('current');
  $(this).addClass('current');
});
