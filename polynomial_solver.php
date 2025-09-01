<?php
// Function to convert from arbitrary base (up to base 36) to decimal number
function baseToDecimal(string $value, int $base): float {
    $digits = '0123456789abcdefghijklmnopqrstuvwxyz';
    $value = strtolower($value);
    $result = 0;

    for ($i = 0; $i < strlen($value); $i++) {
        $digit = strpos($digits, $value[$i]);
        if ($digit === false || $digit >= $base) {
            throw new Exception("Invalid digit '{$value[$i]}' for base {$base}");
        }
        $result = $result * $base + $digit;
    }
    return $result;
}

// Gaussian elimination to solve linear system V * A = Y
function gaussianElimination(array $matrix, array $values): array {
    $n = count($values);

    // Forward elimination
    for ($i = 0; $i < $n; $i++) {
        // Pivot selection
        $maxEl = abs($matrix[$i][$i]);
        $maxRow = $i;
        for ($k = $i + 1; $k < $n; $k++) {
            if (abs($matrix[$k][$i]) > $maxEl) {
                $maxEl = abs($matrix[$k][$i]);
                $maxRow = $k;
            }
        }

        // Swap max row with current row
        if ($maxRow != $i) {
            $tmp = $matrix[$i];
            $matrix[$i] = $matrix[$maxRow];
            $matrix[$maxRow] = $tmp;

            $tmpVal = $values[$i];
            $values[$i] = $values[$maxRow];
            $values[$maxRow] = $tmpVal;
        }

        // Eliminate below
        for ($k = $i + 1; $k < $n; $k++) {
            $c = -$matrix[$k][$i] / $matrix[$i][$i];
            for ($j = $i; $j < $n; $j++) {
                if ($i == $j) {
                    $matrix[$k][$j] = 0;
                } else {
                    $matrix[$k][$j] += $c * $matrix[$i][$j];
                }
            }
            $values[$k] += $c * $values[$i];
        }
    }

    // Back substitution
    $x = array_fill(0, $n, 0);
    for ($i = $n - 1; $i >= 0; $i--) {
        $x[$i] = $values[$i] / $matrix[$i][$i];
        for ($k = $i - 1; $k >= 0; $k--) {
            $values[$k] -= $matrix[$k][$i] * $x[$i];
        }
    }

    return $x;
}

// Read the JSON input file (make sure this file exists)
$json = file_get_contents('testcase.json');
$data = json_decode($json, true);

if ($data === null) {
    die("Invalid JSON input or file not found\n");
}

// Extract keys n and k
$n = intval($data['keys']['n']);
$k = intval($data['keys']['k']);
$m = $k - 1; // degree of polynomial

// Prepare arrays for x and y values
// Assume x = 1, 2, ..., k for roots to build Vandermonde matrix
$x_vals = [];
$y_vals = [];

for ($i = 1; $i <= $k; $i++) {
    if (!isset($data[(string)$i])) {
        die("Missing root $i in input\n");
    }
    $base = intval($data[(string)$i]['base']);
    $val = $data[(string)$i]['value'];
    $decimal = baseToDecimal($val, $base);

    $x_vals[] = $i;
    $y_vals[] = $decimal;
}

// Build Vandermonde matrix (k x k)
$V = [];
for ($i = 0; $i < $k; $i++) {
    $row = [];
    for ($j = 0; $j <= $m; $j++) {
        $row[] = pow($x_vals[$i], $j);
    }
    $V[] = $row;
}

// Solve system for coefficients
$coefficients = gaussianElimination($V, $y_vals);

// Print polynomial coefficients
foreach ($coefficients as $coef) {
    echo $coef . "\n";
}
