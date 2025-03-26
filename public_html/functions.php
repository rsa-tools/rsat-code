<?php
    
    

  // NeAT TITLE
Function title($title) {
  echo "<span class='menu-toggle'><a href='#' id='menu-toggle'><i class='fa fa-bars fa-2x' style='color:green;'></i></a></span><H3><a href='NeAT_home.html'>NeAT</a> - $title</H3>\n";
  }


// NeAT ERROR
Function error($error) {
  echo "<H4>Error</H4>\n<blockquote class='error'>$error</blockquote></h4><br>";
}


// NeAT WARNING
Function warning($warning) {
  echo "<H4>Warning</H4>\n<blockquote class ='warning'>$warning</blockquote><br>";
}


// NeAT WARNING
Function demo($demo) {
  echo "<H4>Comment on the demonstration example</H4>\n<blockquote class ='demo'>$demo</blockquote><hr>\n\n";
}


// NeAT INFO
Function info($info) {
  echo "<h4>Info</h4><class='information'><font color='#00BB00'>".str_replace($properties['RSAT'], '$RSAT', $info)."</font></class><br>\n\n";
  //    echo "<h2>Info</h2><blockquote class='information'><font color='#00BB00'>$info</font></blockquote><br>";
}


// NeAT INFO LINK
// This function adds an hypertext link
Function info_link($info, $link_url) {
  echo "<H4>Info: </H4><blockquote class = 'info'><a href = '$link_url'>$info </a></blockquote></a><br>";
}


Function uploadFile($file) {
  //     $tmpDir = "tmp/";
  //     $tmpFile = $_FILES[$file]["name"];
  //     $now = date("Ymd_His");
  //     $tmpFile = $tmpFile.$now;
  //     $tmpFileName = $tmpDir.$tmpFile;
   $tmpFileName = getTempFileName('upload');
  //   echo ("TEMP $tmpFileName");
  if (is_uploaded_file($_FILES[$file]['tmp_name'])) {
    if (rename($_FILES[$file]['tmp_name'], $tmpFileName)) {
      //             echo "File ".$_FILES[$file]['tmp_name']." was moved to  $tmpDir/$tmpFile <br>";
    } else {
      echo "Could not move $_FILES[$file]['tmp_name']"." check that $tmpDir exists<br>";
    }
  } else {
    echo "File ".$_FILES[$file]['tmp_name']." could not be uploaded<br>";
    //     print_r($_FILES[$file]);
  }
  return $tmpFileName;
}


Function getTempFileName($prefix, $ext) {
  // Main directory
  global $tmp;
  if (isset($tmp)) {
    $tmpDir = $tmp."/";
  } else {
    $tmpDir = "tmp/";
  }

  // Append user name
  global $_ENV; // ceci ne marche pas
#  $user = $_ENV["USER"]; // BUG: ceci ne marche pas
#  $user = getenv('REMOTE_USER'); // ceci ne marche pas
#  $user = getenv('USER'); // ceci ne marche pas
#  $processUser = posix_getpwuid(posix_geteuid()); // ceci ne marche pas sur
               #  RSAT mais bien sur mon portable, je suppose qu''il faut
               #  installer posix, mais je prefere eviter car je devrais le
               #  faire sur tous les miroirs
#  $user = $processUser['name'];
  $user = `whoami`; // THIS WORKS BY IT IS TOO TRICKY TO CLAL A SYSTEM COMMAND.I should check why all the solutions above do not work
  $user = rtrim($user); ##remove space at the end of username ADDED By Didier Croes 11/04/2014
  $tmpDir .= $user."/";

  // Append
  $tmpDir .=  date("Y/m/d/");
  mkdir($tmpDir, 0777, TRUE); // Create temporary subdir if it does not exist
#  $tmpFile = $_FILES[$file]["name"];
#  $tmpFile = $tmpFile.$prefix;
  $now = date("Ymd_His");
  $suffix = randchar(4);  
  #  $suffix = mktemp('XXXXX');
  $tmpFile = $prefix.'_'.$suffix."_".$now;
  if (isset($ext)) {
    $tmpFile .= $ext;
  }
  //  print_r("<br>tmpDir = ".$tmpDir."</br>");
  //  print_r("<br>tmpFile = ".$tmpFile."</br>");

  // return the complete path to the temporary file
  return $tmpDir.$tmpFile;
}




Function storeFile($file) {
  global $properties;
  # transforms the $RSAT variable to the real RSAT path
  if (isset($properties['RSAT'])) {
    $file = str_replace("\$RSAT", $properties['RSAT'], $file);
  }
  $fh = fopen($file, 'r');
  $theData = "";
  if (filesize($file) > 0) {
    $theData = fread($fh, filesize($file));
    fclose($fh);
  }
  return $theData;
}


Function writeTempFile($prefix, $data) {
  $tmpFileName = getTempFileName($prefix);

#  $temp_command = "mktemp tmp/$prefix.XXXXX";
  $temp_command = "mktemp ".$tmpFileName;
  exec($temp_command, $temp_file);
  $temp_file = $temp_file[0];
  $temp_file = trim_text($temp_file);
  $temp_file = trim($temp_file);
  $fh = fopen($temp_file, 'w');
  fwrite($fh, $data);
  fclose($fh);
  return $temp_file;
}


/**
 * Strip special invisible characters
 */
Function trim_text($text) {
  $trimmed_text = "";
  $lines = explode("\n",$text);
  $array_count = count($lines);
  for($y=0; $y<$array_count; $y++) {
    $trimmed_text .= rtrim($lines[$y])."\n";
  }
  return $trimmed_text;
}


/**
 * Echo a command by suppress the full path to RSAT
 */
Function store_command($command, $name, $cmd_handle) {
  $clean_command = preg_replace('/(\')*\S+rsa\-tools/', '\\1\$RSAT', $command);
  $clean_command = preg_replace('/\/+/', '/', $clean_command);
  $clean_command = preg_replace('/\'(\S+)\'/', '\\1', $clean_command);
  //  echo ("<p><b>$name:</b> $clean_command</p>");
  fwrite($cmd_handle, "# ".$name."\n");
  fwrite($cmd_handle, $clean_command."\n\n");
}


/**
 * Print a table with URLs to the result files
 */
Function print_url_table($URL) {
  echo "<table class=\"resultlink\">\n";
  echo "<tr><th colspan='2'>Result file(s)</th></tr>\n";
  foreach ($URL as $key => $value) {
    echo "<tr><td>",$key,"</td><td><a href = '",$value,"'>",$value,"</a></td></tr>\n"; 
  }
  echo "</table>\n";
  echo"<hr>\n";
}


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
    if (!preg_match("/^\#/", $line)) {
      $property = explode('=', $line);
      if (isset($property[1])) {
	$prop_array[$property[0]] = $property[1];
      }
    }
  }
  return $prop_array;
}

// Requires parsing file $RSAT/version.txt, returns 1st word in file $RSAT/version.txt,
// which should be an ISO data such as 2024-08-28
Function getGitLastReleaseDate(){

   // Guess the main RSAT directory
   $rsat_main = getcwd()."/..";

   $date = shell_exec("head $rsat_main/version.txt");

   return chop($date);
}

// Deprecated for the problems described and the fact that Docker and conda versions
// are built from releases not repos, see https://github.com/rsa-tools/rsat-code/issues/49
//
// PHP Warning:  Trying to access array offset on value of type null in 
// /var/www/html/rsat/public_html/functions.php on line ...
// fatal: detected dubious ownership in repository at '/var/www/html/rsat'
// To add an exception for this directory, call:
//   git config --global --add safe.directory /var/www/html/rsat
//
// explanation: user www-data does not own repository and cannot do git pull,
//              observed in Debian 12 with php 8.2
//
Function getGitLastCommitDate(){
	$date = shell_exec('git log | head -n 4 | grep Date'); 
	//Date:   Thu Jun 15 10:52:32 2023 +0200

	if (isset($date)) {
		$date = preg_replace('/^\s*\w+\s+/', '', $date);
		$date = preg_replace('/\D+\w+\s*$/', '', $date);

	} else {
		$rsat_path = getcwd() . "/../";
		$date = date("F d H:i:s Y e", filemtime( $rsat_path . '.git'));
	}

	return $date;
}

////////////////////////////////////////////////////////////////
// Operations done when loading each php page

// Guess the main RSAT directory
$rsat_main = getcwd()."/..";
    $x = getcwd();
    if(strpos($x, "/tutorials") !== false){
        $x = str_replace("/tutorials", "", $x);
        $rsat_main = $x."/..";
    }
// log file
$rsat_logs = $rsat_main."/logs";

// Load properties from the site-specific configuration file RSAT_config.props
$properties = load_props($rsat_main."/RSAT_config.props");
$properties['git_release'] = getGitLastReleaseDate();
//$properties['git_date'] = getGitLastCommitDate();//deprecated

// Check if the property rsat_www has been set to "auto", and, if so,
// set it automatically
if ($properties['rsat_www'] == "auto") {
  $ip = $_SERVER['HTTP_HOST'];
  // print_r('IP = '.$ip."<p>\n");
  $properties['rsat_www'] =  'http://'.$ip.'/rsat/';
}


// Check if the property neat_www_root is set and if it is set to "auto"
if (!isset($properties['neat_www_root']) || $properties['neat_www_root'] == "auto") {
  $ip = $_SERVER['HTTP_HOST'];
  //print_r('IP = '.$ip);
  $properties['neat_www_root'] =  'http://'.$ip.'/rsat/';
}

@ $neat_www_root = $properties['neat_www_root'];

// Some particular properties are defined as global variables, for convenience
$tmp = $properties['rsat_tmp'];
$RSAT = $properties['RSAT'];
$rsat_www = $properties['rsat_www'];
$log_name = $properties['rsat_site'];
date_default_timezone_set("Europe/Paris");

// automatic definition of property neat_www_root
if ($properties['neat_www_root'] == "auto") {
  $properties['neat_www_root'] = $rsat_www;
}
$neat_www_root = $properties['neat_www_root'];

// automatic definition of neat_ws
if (!isset($properties['neat_ws']) || $properties['neat_ws'] == "auto") {
//if ($properties['neat_ws'] == "auto") {
  $properties['neat_ws'] = $rsat_www;
}

// address of the WSDL document
$neat_wsdl = $properties['neat_ws']."/web_services/RSATWS.wsdl";

// // Karoline: property neat_java_host not required for me, just need the host name here
// // in future: obtain host from url address
// // $scheme = parse_url($rsat_www,PHP_URL_SCHEME);
// // $host = parse_url($rsat_www,PHP_URL_HOST);
// // $neat_java_host = $scheme."://".$host;
// // for the moment: java tools run only on ulb host
// $neat_java_host = $properties['neat_java_host'];
// // host may include tomcat port
// if (isset($properties['tomcat_port'])) {
//   $tomcat_port = $properties['tomcat_port'];
//   if (strcmp($tomcat_port,"") != 0){
//     $neat_java_host = $neat_java_host . ":" . $tomcat_port;
//   }
//   $neat_java_wsdl = $neat_java_host . "/be.ac.ulb.bigre.graphtools.server/wsdl/GraphAlgorithms.wsdl";
//   $neat_java_remote_wsdl = $neat_java_host . "/be.ac.ulb.bigre.graphtools.server/wsdl/GraphAlgorithmsNeAT.wsdl";
// }

// LOG
$year = date("Y");
$month = date("m");
$neat_log_file = sprintf ("$rsat_logs/log-file_$log_name"."_neat_%04d_%02d", $year,$month);
$rsat_log_file = sprintf ("$rsat_logs/log-file_$log_name"."_%04d_%02d", $year,$month);


// This function returns the name of the script executing it
Function getCurrentScriptName() {
  $currentFile = $_SERVER["SCRIPT_NAME"];
  $parts = Explode('/', $currentFile);
  $currentFile = $parts[count($parts) - 1];
  return $currentFile;
}


// This function replaces all spaces of a string by tabulation
// If the line starts with a ';' or a '#' it is skipped.
Function space_to_tab($string) {
  $result = "";
  $lines = explode("\n",$string);
  $array_count = count($lines);
  for($i=0; $i<$array_count; $i++) {
    $line = $lines[$i];
    if (!preg_match("/^\#/", $line) && !preg_match("/^\;/", $line)) {
      $line_sp = str_replace(" ", "\t", $line);
      $result .= $line_sp."\n";
    } else {
      $result .= $line."\n";
    }
  }
  return $result;
}


// This function replaces the given number of spaces in a row inside a string by tabulation
// If the line starts with a ';' or a '#' it is skipped.
// This function allows to handle input where identifiers may contain less than the given
# number of spaces in a row.
Function spaces_to_tab($string, $num) {
  $result = "";
  $replace = "";
  for($i=0; $i<$num; $i++){
    $replace = $replace." ";
  }
  $lines = explode("\n",$string);
  $array_count = count($lines);
  for($i=0; $i<$array_count; $i++) {
    $line = $lines[$i];
    if (!preg_match("/^\#/", $line) && !preg_match("/^\;/", $line)) {
      $line_sp = str_replace($replace, "\t", $line);
      $result .= $line_sp."\n";
    } else {
      $result .= $line."\n";
    }
  }
  return $result;
}


# This function converts a file name from its complete path
# its URL on the RSAT webserver
# For example: /home/rsat/rsa-tool/public_html/tmp/brol.truc
# will be converted to a link to the corresponding URL on the web server.
Function rsat_path_to_url ($file_name) {
  global $rsat_www;
  global $properties;
  # transforms the $RSAT variable to the real RSAT path
  $file_name = str_replace("\$RSAT", $properties['RSAT'], $file_name);


#  $temp_file = rtrim($file_name);
#  $temp_file = explode('/',$temp_file);
#  $temp_file = end($temp_file);
#  $resultURL = $rsat_www."/tmp/".$temp_file;
  $resultURL = str_replace($properties['RSAT']."/public_html", $properties['rsat_www'], $file_name);
  return $resultURL;
}


// This function returns the name of the script executing it
Function AlphaDate() {
  $my_date = exec("date +%Y_%m_%d.%H%M%S", $my_date);
  trim($my_date);
  return $my_date;
}


// This function returns the name of the script executing it
Function check_integer($string) {
  return (preg_match("/[0-9]*/", $string));
}


////////////////////////////////////////////////////////////////
/* 

Store info into a log file in a conveninent way for subsequent login
statistics.

Usage:
     UpdateLogFile();
     UpdateLogFile($script_name);
     UpdateLogFile($script_name, $message);

If script name is empty or null, the program determines the name of
the file.
*/
Function UpdateLogFile($suite ,$script_name, $message) {
  if ($script_name == "") {
    $script_name = getCurrentScriptName();
  }

  // LOAD GLOBAL VARIABLES
  global $log_file, $log_name, $rsat_log_file, $neat_log_file;
  if ($suite == "rsat") {
    $log_file = $rsat_log_file;
  } else {
    $log_file = $neat_log_file;
  }
  // LOAD OTHER VARIABLES
  $my_date = AlphaDate();
  $user = getenv('REMOTE_USER');
  $address = getenv('REMOTE_ADDR');
  $host = getenv('REMOTE_HOST');
  $e_mail = "";
  $user_address_at_host = $user."@".$address." (".$host.")";
  $to_write = $my_date."\t".$log_name."\t".$user_address_at_host."\t".$script_name."\t".$e_mail."\t".$message."\n";

  // Write to the file
  $log_handle = fopen($log_file, 'a');
  fwrite($log_handle, $to_write);
  fclose($log_handle);
  chmod ($log_file, 0777);
}


Function checkNeatTutorial($tutorial_url) {
  global $rsat_www;
  $pdf_tutorial = $rsat_www.'/tutorials/neat_tutorial.pdf';
  $address = $pdf_tutorial;
  if (file_exists($tutorial_url)) {
    $address = $tutorial_url;
  }
  return $address;
}


Function HourglasOn() {
  echo("<div id='hourglass' class='hourglass'><img src='images/animated_hourglass.gif' height='50' border='1'></div>");
}


Function hourglass($status) {
  if ($status == "on") {
    echo("<div id='hourglass' class='hourglass'><img src='images/animated_hourglass.gif' height='50' border='1'></div>");
  } else {
    echo("<div id='hide' class='hide'><img src='images/hide_hourglass.jpg' height='60' border='0'></div>");
  }
}


function urlfilesize($url,$thereturn) {
  if (substr($url,0,4)=='http') {
    $x = array_change_key_case(get_headers($url, 1),CASE_LOWER);
    $x = $x['content-length'];
  }
  else {
    $x = @filesize($url); }
  if (!$thereturn) {
    return $x ;
  }
  elseif($thereturn == 'mb') {
    return round($x / (1024*1024),2) ;
  }
  elseif($thereturn == 'kb') {
    return round($x / (1024),2) ;
  }
}


function randchar($length) {
  $result = "";
  $abc= array("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9");
  for ($i = 0; $i < $length; $i++) {
    $result .= $abc[rand(0,(count($abc)-1))];
  }
  return $result;
}

function converttab2table($tabdelimtable, $tableprops) {
	//class=sortable id=seedstable
	$output = "<table ".$tableprops." >\n";
	$linearray = preg_split ("#\n#",$tabdelimtable);
	foreach ($linearray as $line){
		if(startsWith($line,"#",false)){
			$line ="<thead><tr><th> </th><th>".$line;
			$line = str_replace("\t","</th><th>",$line );
			$line .="</th></tr></thead>\n";
		}else{
			$contentarray = preg_split ("#\t#",$line);
			$line ="<tr><td><input name=chkID id=chkID type=checkbox  value=\"".$contentarray[0]."\" checked=checked></td><td>".$line;
			$line = str_replace("\t","</td><td>",$line );
			$line .="</td></tr>\n";
		}
		$output.=$line; 
	}
	$output.="</table>\n";
	return $output;
}

function startsWith($haystack,$needle,$case=true)
{
   if($case)
       return strpos($haystack, $needle, 0) === 0;

   return stripos($haystack, $needle, 0) === 0;
}

function endsWith($haystack,$needle,$case=true)
{
  $expectedPosition = strlen($haystack) - strlen($needle);

  if($case)
      return strrpos($haystack, $needle, 0) === $expectedPosition;

  return strripos($haystack, $needle, 0) === $expectedPosition;
}

Function printMenu($menu){
        if(strcmp($menu, "RSAT") == 0){
            ob_start(); // begin collecting output
            include 'menu.php';
            $result = ob_get_clean(); 
            echo($result);
        }else{
            ob_start(); // begin collecting output
            include 'menu_graph.php';
            $result = ob_get_clean(); 
            echo($result);			
        }
    }

?> 
<?php
ini_set('max_execution_time', 2400);
ini_set('max_input_time', 6000);
ini_set('memory_limit', "100M");
ini_set('display_errors', 'on');
ini_set('post_max_size', 80000000);
ini_set('error_reporting', "E_ALL & ~E_NOTICE");
ini_set('upload_max_filesize', 100000000);
?>
