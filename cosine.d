/* Calculates a cosine score between two arrays of equal length.
 *
 * Author: Jonathan Samson
 * Date: 03 August 2020
 */
module cosine;
import std.math;
import std.exception;

private real find_numerator(real[] a, real[] b)
/* Calculates the numerator for the cosine value.
 * Arguments:
 *	a - The first vector to be compared.
 *	b - The second vector to be compared.
 * Returns:
 *	sum - The overall sum of AiBi.
 */
{
	real sum = 0.0;
	for(int i = 0; i < a.length; ++i)
	{
		sum += (a[i] * b[i]);	
	}
	return sum;
}
unittest
{
	real[] a = [1, 1, 1, 1];
	real[] b = [1, 1, 1, 1];
	assert(find_numerator(a, b) == 4);
}

private real find_denominator(real[] a, real[] b)
/* Calculates the denominator for the cosine value.
 * Arguments:
 *	a - The first vector to be compared.
 *	b - The second vector to be compared.
 * Returns:
 *	product - The overall denominator of the consine value.
 */
{
	real sum_a = 0.0;
	real sum_b = 0.0;
	foreach(value; a)
	{
		sum_a += value ^^ 2;
	}
	foreach(value; b)
	{
		sum_b += value ^^ 2;
	}
	real product = sqrt(sum_a) * sqrt(sum_b);
	return product; 
}
unittest
{
	real[] a = [3, -4];
	real[] b = [5, 12];
	assert(find_denominator(a, b) == 65);
}

real calculate_cosine_score(real[] a, real[] b)
/* Calculates the cosine score between the two arrays of equal length.
 * Arguments:
 *	a - The first vector to be compared.
 *	b - The second vector to be compared.
 * Returns:
 *	cosine - The cosine score between a and b.
 */
{
	enforce(a.length == b.length, "The arrays must be of equal length.");
	real cosine = find_numerator(a, b) / find_denominator(a, b);
	return cosine;
}
unittest
{
	real[] a = [3, -4];
	real[] b = [3, -4];
	assert(calculate_cosine_score(a, b) == 1);

	b = [-3, 4];
	assert(calculate_cosine_score(a, b) == -1);
	
	b = [4, 3];
	assert(calculate_cosine_score(a, b) == 0);

	assert(calculate_cosine_score(a, b) == 
	       calculate_cosine_score(b, a));

	b = [1, 2, 3];
	assertThrown(calculate_cosine_score(a, b));
}

void main()
{
	return;
}
