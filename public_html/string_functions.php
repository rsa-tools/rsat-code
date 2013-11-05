<?php
// List of functions dedicated to the treatment of STRING database via REST webservices using NeAT.
?>
<?php
  // Check if the gene names exists in STRING
  // return their description in an array
//   // 1 String ID
  // 2 Usual name
  // 3 annotation
  Function resolveName($names, $organism_id) {
    $result = array ();
    for ($i = 0; $i < count($names); $i++) {
      $name = $names[$i];
      if ($name == "") {
        continue;
      }
      $tempfile = writeTempFile("resolveString", "");
      $wget_command = "wget 'http://stitch.embl.de/api/tsv-no-header/resolve?identifier={$name}&species=$organism_id&echo_query=1' -O $tempfile";
      
//       $wget_command = "wget 'http://string80.embl.de/api/tsv/interactors?identifier=$name&species=$organism_id' -O $tempfile";
      exec ($wget_command);
      $resolve_content = storeFile ($tempfile);
      $lines = explode("\n", $resolve_content);
      $result_size = count ($result);
      for ($j = 0; $j < count ($lines); $j++) {
        $line = $lines[$j];
        $linecp = explode("\t", $line);
        if ($line == "") {
          continue;
        }
        if ($linecp[0] == 'not found') {
          $result[$result_size][1] = "$name";
          $result[$result_size][0] = "Not found in STRING";
          $result[$result_size][2] = "No available definition";
          break;
        }
        $result[$result_size][0] = $linecp[4]; 
        $result[$result_size][1] = $linecp[1];
        $result[$result_size][2] = $linecp[5];
        $result_size++;
      }
    }
    return $result;
  }
?>

<?php
  // Returns an array containing
  // 1) the organism ID
  // 2) the usual name of the organism
  Function readStringOrganisms() {
    $organisms = file_get_contents("http://string.embl.de/newstring_download/species.v9.0.txt");
    $lines = explode("\n",$organisms);
    $sorted_lines = array();
    $array_count = count($lines);
    for($y=0; $y<$array_count; $y++) {
      $line = trim($lines[$y]);
      if (!preg_match("/^\#\#/", $line)) {
        $linecp = explode("\t", $line);
        array_push($sorted_lines, $linecp[2]."\t".$linecp[0]);
      }
    }
    sort($sorted_lines);

    $array_count = count($sorted_lines);
    for($y=0; $y<$array_count; $y++) {
      $line = trim($sorted_lines[$y]);
      
      if (!preg_match("/^\#\#/", $line)) {
        $linecp = explode("\t", $line);
//         print_r($linecp);
        $organism_array[$y][0] = $linecp[1];
        $organism_array[$y][1] = $linecp[0];
      }
    }
    return $organism_array;
  }
?>