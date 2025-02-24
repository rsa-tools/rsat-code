<html xmlns="http://www.w3.org/1999/xhtml">
	<head>
<?php

if ($properties["GOOGLE_ID"]){
	$google_id = $properties["GOOGLE_ID"];
	echo '<!-- Google tag (gtag.js) -->';
	echo '<script async src="https://www.googletagmanager.com/gtag/js?id=$google_id"></script>';
	echo '<script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag("js", new Date())';
	echo "gtag('config', '$google_id');"
	echo '</script>';
}

if ($result) {
	$res = " - results";
	$class = "results";
	echo '<link rel="stylesheet" type="text/css" href = "css/main_grat.css" media="screen">',"\n";
} else  {
	$res = "";
	$class = "form";
	echo '<link rel="stylesheet" type="text/css" href="css/main.css" media="screen,projection,print"/>',"\n";
}
?>
<!--  	<link rel="stylesheet" type="text/css" href="tabs.css" media="screen,projection,print"/>-->	
		<script src="js/RSAT_menu.js" type="text/javascript"></script>
		<script src="js/RSAT_tabs.js" type="text/javascript"></script>
		<title>RSAT - <?php echo $prog_name, $res; ?></title>
		<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	</head>

	<body class="<?php echo $class; ?>">
		<h3 align='center'><a href="<?php echo $properties['rsat_www']?>">RSAT</a> - <?php echo $prog_name, $res; ?></h3><br/>
