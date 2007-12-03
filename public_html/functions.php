<?php
  // Grat TITLE
  Function title($title) {
    echo "<H3><a href='NeAT_home.html'>NeA-tools</a> - $title</H3>\n";
  } 
?>

<?php
  // Grat ERROR
  Function error($error) {
    echo "<H4>Error : </H4><blockquote class='error'>$error</blockquote></h4><br>";
  } 
?>

<?php
  // Grat WARNING
  Function warning($warning) {
    echo "<H4>Warning : </H4><blockquote class ='warning'>$warning</blockquote><br>";
  } 
?>


<?php
  // Grat WARNING
  Function demo($demo) {
    echo "<H4>Remark : </H4><blockquote class ='demo'>$demo</blockquote><br>";
  } 
?>

<?php
  // Grat INFO
  Function info($info) {
    echo "<H4>Info : </H4><blockquote class = 'info'>$info </blockquote><br>";
  } 
?>

<?php
    Function uploadFile($file) {
    $repertoireDestination = "tmp/";
    $nomDestination = $_FILES[$file]["name"];
    $now = date("Ymd_His");
    $nomDestination = $nomDestination.$now;
    
    if (is_uploaded_file($_FILES[$file]['tmp_name'])) {
        if (rename($_FILES[$file]['tmp_name'], $repertoireDestination.$nomDestination)) {
//             echo "File ".$_FILES[$file]['tmp_name']." was moved to  $repertoireDestination/$nomDestination <br>";
        } else {
            echo "Could not move $_FILES[$file]['tmp_name']"." check that $repertoireDestination exists<br>";
        }          
    } else {
       echo "File ".$_FILES[$file]['tmp_name']." could not be uploaded<br>";
    }
    return $repertoireDestination.$nomDestination;
  }
?>

<?php
    Function storeFile($file) {
      $fh = fopen($file, 'r');
      $theData = fread($fh, filesize($file));
      fclose($fh);
      return $theData;
  }
?>

<?php
/**
 * Strip special invisible characters
 */
Function trim_text($text) {
   $trimmed_text = "";
   $lines = explode("\n",$text);
   $array_count = count($lines);
   for($y=0; $y<$array_count; $y++) {
     $trimmed_text .= trim($lines[$y])."\n";
   }
   return $trimmed_text;
}
?>
<?php
/**
 * Read a property file $props and return a hash
 */
Function load_props($props) {
  $prop_array = array();
  $properties = storeFile($props);
  $lines = explode("\n",$properties);
  $array_count = count($lines);
  for($y=0; $y<$array_count; $y++) {
   $line = trim($lines[$y]);
   $property = explode('=', $line);
   $prop_array[$property[0]] = $property[1];
  }
  return $prop_array;
}
?>


<?php
  # LOAD PROPERTIES
  $rsat_public_html = getcwd()."/..";;
  $properties = load_props($rsat_public_html."/RSAT_config.props");
  $tmp = $properties[rsa_tmp];
  $WWW_RSA = $properties[www_rsa];
//   ini_set('soap.wsdl_cache_enabled', 0);
?>
