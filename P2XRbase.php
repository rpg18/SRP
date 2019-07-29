<?php

// connect
$connection = db_connect();

//define gene
$target = mysqli_real_escape_string($connection, $_POST['gene']);

// define SQL
$sql = "SELECT * FROM website WHERE Gene = '".$target."'";

// Attempt select query execution
if($result = mysqli_query($connection, $sql)){
   if(mysqli_num_rows($result) > 0){
      echo "<table>";
         echo "<tr>";
            echo "<th>Gene</th>";
            echo "<th>Run</th>";
            echo "<th>age</th>";
            echo "<th>cell_type</th>";
            echo "<th>tissue</th>";
            echo "<th>k9_cluster</th>";
            echo "<th>Gene_Count</th>"; 
         echo "</tr>";
      while($row = mysqli_fetch_array($result)){
         echo "<tr>";
            echo "<td>" . $row['Gene'] . "</td>";
            echo "<td>" . $row['Run'] . "</td>";
            echo "<td>" . $row['age'] . "</td>";
            echo "<td>" . $row['cell_type'] . "</td>";
            echo "<td>" . $row['tissue'] . "</td>";
            echo "<td>" . $row['k9_cluster'] . "</td>";
            echo "<td>" . $row['Gene_Count'] . "</td>";
         echo "</tr>";
      }
      echo "</table>";
      // Free result set
      mysqli_free_result($result);
   } 
   else{
        echo "No records matching your query were found.";
   }
} 
else{
    echo "ERROR: Could not execute $sql. " . mysqli_error($connection);
}



function db_connect() {
   //static variable, avoid multiple connections
   static $connection;
	
   // connect unless connection exists
   if(!isset($connection)) {
      // load config data
      $config = parse_ini_file('p2xdb.ini');	
      // Try and connect to the database
      $connection = mysqli_connect('localhost',$config['username'],$config['password'],$config['dbname']);
   }
   // If connection was not successful, handle the error
   if($connection === false) {
        echo "connection fail";
      // Handle error 
      return mysqli_connect_error();	
   }
   return $connection;
}

?>
