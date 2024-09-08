from collections.abc import Sequence
from itertools import chain

def separate_string(input_string: str) -> list[str]:
    """
    Separates the input string by newline characters or semicolons and returns a list of the substrings.

    Parameters
    ----
    input_string (str): The string to be separated.

    Returns
    -------
    list: A list of substrings separated by newline characters or semicolons.

    Examples
    -------
    >>> input_string = "Hello world;\nThis is a test;\nAnother line"
    >>> print(separate_string(input_string))
    Output: ['Hello world', 'This is a test', 'Another line']
    """
    # Splitting the string by either newline characters or semicolons
    lines = input_string.splitlines()  # Split by newlines first

    # Further split each line by semicolons
    substrings = chain.from_iterable([line.split(';') for line in lines])

    # Filter out any empty strings that may result from splitting
    return [substring.strip() for substring in substrings if substring]

def remove_comments(input_string: str) -> str:
    """
    Removes comments from the input string.

    Args:
    input_string (str): The string from which comments need to be removed.

    Returns:
    str: A string with comments removed, where comments start with a '#' and end at the next newline or semicolon.

    Examples:
    >>> input_string = "1+1#sum\n#new var\na=1"
    >>> print(remove_comments(input_string))
    Output: "1+1\n\na=1"
    """
    result = []
    i = 0
    while i < len(input_string):
        if input_string[i] == '#':
            # Skip characters until a newline or semicolon is found
            while i < len(input_string) and input_string[i] not in ('\n', ';'):
                i += 1
            # If it's a semicolon, skip it as well
            if i < len(input_string) and input_string[i] == ';':
                i += 1
        if i < len(input_string):
            result.append(input_string[i])
        i += 1
    
    return ''.join(result)

def simplify_terms(terms: Sequence[tuple[int | float, str]]) -> list[tuple[int | float, str]]:
    """
    Simplify the reactions terms.

    Parameters
    ----------
    terms : Sequence[tuple(int | float, str)]
        The terms to simplify.

    Returns
    -------
    list[tuple(int | float, str)]
        The simplified terms.

    Details
    -------
    The simplification is made in the following way:
        - If the coefficient is absent, it is set to 1.
        - If the coefficient is 0, the term is removed.
        - If two or more terms have the same species, the coefficient is summed.
    """
    term_dict = {}
    
    for coefficient, species in terms:
        if species in term_dict:
            term_dict[species] += coefficient
        else:
            term_dict[species] = coefficient
    
    # Remove any terms with a coefficient of 0 and convert to list of tuples
    simplified_terms = [(coef, spec) for spec, coef in term_dict.items() if coef != 0]
    
    return simplified_terms

def str2int_or_float(string: str) -> int | float:
    try:
        return int(string)
    except ValueError:
        return float(string)