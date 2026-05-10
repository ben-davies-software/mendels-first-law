"""
Mendel's First Law Solver
=========================
Rosalind Problem: IPRB

Given a population of k homozygous dominant (AA), m heterozygous (Aa),
and n homozygous recessive (aa) organisms, calculate the probability that
two randomly selected mating organisms will produce an offspring displaying
the dominant phenotype.

Author: Ben Davies 
Created for Rosalind bioinformatics platform
"""

from fractions import Fraction
from typing import Union


class Population:
    """Represents a population of organisms with Mendelian genotypes."""

    GENOTYPES = ("homozygous dominant (AA)", "heterozygous (Aa)", "homozygous recessive (aa)")

    def __init__(self, k: int, m: int, n: int) -> None:
        """
        Initialise the population.

        Args:
            k: Number of homozygous dominant individuals (AA)
            m: Number of heterozygous individuals (Aa)
            n: Number of homozygous recessive individuals (aa)

        Raises:
            ValueError: If any value is negative or total population is fewer than 2.
        """
        self._validate(k, m, n)
        self.k = k  # AA
        self.m = m  # Aa
        self.n = n  # aa

    @staticmethod
    def _validate(k: int, m: int, n: int) -> None:
        """Validate population counts."""
        for label, value in [("k", k), ("m", m), ("n", n)]:
            if not isinstance(value, int):
                raise TypeError(f"Value for '{label}' must be an integer, got {type(value).__name__}.")
            if value < 0:
                raise ValueError(f"Value for '{label}' must be a non-negative integer, got {value}.")

        total = k + m + n
        if total < 2:
            raise ValueError(
                f"Total population must be at least 2 to select a mating pair, "
                f"but got {total} (k={k}, m={m}, n={n})."
            )

    @property
    def total(self) -> int:
        """Total number of organisms in the population."""
        return self.k + self.m + self.n

    def __repr__(self) -> str:
        return (
            f"Population(k={self.k} [AA], m={self.m} [Aa], n={self.n} [aa], "
            f"total={self.total})"
        )

    def summary(self) -> str:
        """Return a human-readable summary of the population."""
        lines = [
            "Population Summary",
            "==================",
            f"  Homozygous dominant (AA) : {self.k}",
            f"  Heterozygous        (Aa) : {self.m}",
            f"  Homozygous recessive(aa) : {self.n}",
            f"  Total                    : {self.total}",
        ]
        return "\n".join(lines)


class MendelSolver:
    """
    Solves Mendel's First Law (Rosalind IPRB).

    Calculates the probability that two randomly selected organisms from
    a population will produce an offspring with at least one dominant allele
    (dominant phenotype), assuming any two organisms can mate.
    """

    # Offspring dominant probability for each parent pair combination
    # Rows/cols: AA (k), Aa (m), aa (n)
    # Probability that a cross produces a dominant-phenotype offspring:
    #   AA x AA = 1.0   (all offspring AA)
    #   AA x Aa = 1.0   (1/2 AA, 1/2 Aa -- all dominant)
    #   AA x aa = 1.0   (all Aa -- all dominant)
    #   Aa x Aa = 0.75  (1/4 AA, 1/2 Aa, 1/4 aa)
    #   Aa x aa = 0.5   (1/2 Aa, 1/2 aa)
    #   aa x aa = 0.0   (all aa)
    CROSS_DOMINANT_PROB = {
        ("AA", "AA"): Fraction(1, 1),
        ("AA", "Aa"): Fraction(1, 1),
        ("AA", "aa"): Fraction(1, 1),
        ("Aa", "AA"): Fraction(1, 1),
        ("Aa", "Aa"): Fraction(3, 4),
        ("Aa", "aa"): Fraction(1, 2),
        ("aa", "AA"): Fraction(1, 1),
        ("aa", "Aa"): Fraction(1, 2),
        ("aa", "aa"): Fraction(0, 1),
    }

    def __init__(self, population: Population) -> None:
        """
        Initialise solver with a Population instance.

        Args:
            population: A validated Population object.
        """
        if not isinstance(population, Population):
            raise TypeError(f"Expected a Population instance, got {type(population).__name__}.")
        self.population = population

    def _pair_probability(self, count_a: int, count_b: int, same_type: bool) -> Fraction:
        """
        Calculate the probability of selecting a specific ordered pair of organisms.

        For pairs from the same genotype group, the second draw is from the
        remaining pool (without replacement).

        Args:
            count_a: Number of organisms of type A.
            count_b: Number of organisms of type B.
            same_type: True if both organisms are from the same genotype group.

        Returns:
            Probability of selecting this ordered pair.
        """
        total = self.population.total
        if same_type:
            if count_a < 2:
                return Fraction(0)
            return Fraction(count_a, total) * Fraction(count_a - 1, total - 1)
        else:
            return Fraction(count_a, total) * Fraction(count_b, total - 1)

    def calculate(self) -> float:
        """
        Calculate the probability of producing a dominant-phenotype offspring.

        Uses exact arithmetic (Python Fraction) to accumulate the probability
        over all ordered parent-pair combinations, then converts to float.

        Returns:
            Probability as a float, rounded to 5 decimal places.
        """
        k, m, n = self.population.k, self.population.m, self.population.n

        # All ordered pair combinations (with replacement of genotype labels,
        # but without replacement of individual organisms)
        genotype_counts = {"AA": k, "Aa": m, "aa": n}
        genotypes = ["AA", "Aa", "aa"]

        total_prob = Fraction(0)

        for g1 in genotypes:
            for g2 in genotypes:
                same = g1 == g2
                p_pair = self._pair_probability(genotype_counts[g1], genotype_counts[g2], same)
                p_dominant = self.CROSS_DOMINANT_PROB[(g1, g2)]
                total_prob += p_pair * p_dominant

        return float(total_prob)

    def detailed_breakdown(self) -> str:
        """
        Return a detailed breakdown of each cross combination and its contribution.

        Returns:
            Formatted string showing all crosses, their probabilities, and contributions.
        """
        k, m, n = self.population.k, self.population.m, self.population.n
        genotype_counts = {"AA": k, "Aa": m, "aa": n}
        genotypes = ["AA", "Aa", "aa"]

        lines = [
            "",
            "Detailed Cross Breakdown",
            "========================",
            f"{'Cross':<12} {'P(pair)':<18} {'P(dominant|cross)':<22} {'Contribution':<18}",
            "-" * 70,
        ]

        total_prob = Fraction(0)

        for g1 in genotypes:
            for g2 in genotypes:
                same = g1 == g2
                p_pair = self._pair_probability(genotype_counts[g1], genotype_counts[g2], same)
                p_dominant = self.CROSS_DOMINANT_PROB[(g1, g2)]
                contribution = p_pair * p_dominant
                total_prob += contribution

                cross_label = f"{g1} x {g2}"
                lines.append(
                    f"{cross_label:<12} {str(p_pair):<18} {str(p_dominant):<22} {str(contribution):<18}"
                )

        lines.append("-" * 70)
        lines.append(f"{'TOTAL':<12} {'':18} {'':22} {float(total_prob):.5f}")
        return "\n".join(lines)


class InputHandler:
    """Handles user input with validation for the Mendel solver."""

    @staticmethod
    def get_non_negative_int(prompt: str, label: str) -> int:
        """
        Prompt the user for a non-negative integer, with retry on invalid input.

        Args:
            prompt: The prompt string displayed to the user.
            label: The variable label used in error messages.

        Returns:
            A validated non-negative integer.
        """
        while True:
            raw = input(prompt).strip()

            if not raw:
                print(f"  [Error] Input for '{label}' cannot be empty. Please enter a whole number.")
                continue

            try:
                value = int(raw)
            except ValueError:
                print(
                    f"  [Error] '{raw}' is not a valid integer for '{label}'. "
                    f"Please enter a whole number (e.g. 0, 2, 10)."
                )
                continue

            if value < 0:
                print(f"  [Error] '{label}' must be 0 or greater, got {value}.")
                continue

            return value

    @staticmethod
    def get_population_from_user() -> Population:
        """
        Interactively prompt the user for k, m, and n values and return a Population.

        Handles all validation, including ensuring the total population is at
        least 2 for a valid mating pair to exist.

        Returns:
            A validated Population instance.
        """
        print("\nEnter population counts (non-negative integers):")

        while True:
            k = InputHandler.get_non_negative_int("  k (homozygous dominant, AA): ", "k")
            m = InputHandler.get_non_negative_int("  m (heterozygous, Aa)       : ", "m")
            n = InputHandler.get_non_negative_int("  n (homozygous recessive, aa): ", "n")

            try:
                pop = Population(k, m, n)
                return pop
            except ValueError as exc:
                print(f"  [Error] {exc}")
                print("  Please re-enter the values.\n")

    @staticmethod
    def get_population_from_string(data: str) -> Population:
        """
        Parse a space-separated string of three integers into a Population.

        Args:
            data: A string containing three non-negative integers (e.g. "2 2 2").

        Returns:
            A validated Population instance.

        Raises:
            ValueError: If the string cannot be parsed or values are invalid.
        """
        parts = data.strip().split()
        if len(parts) != 3:
            raise ValueError(
                f"Expected exactly 3 space-separated integers (k m n), got {len(parts)}: '{data}'."
            )

        parsed = []
        for label, part in zip(("k", "m", "n"), parts):
            try:
                value = int(part)
            except ValueError:
                raise ValueError(f"Cannot parse '{part}' as an integer for '{label}'.")
            if value < 0:
                raise ValueError(f"Value for '{label}' must be non-negative, got {value}.")
            parsed.append(value)

        return Population(*parsed)


def run_interactive() -> None:
    """Run the solver in interactive mode, prompting the user for input."""
    print("=" * 60)
    print("  Mendel's First Law Solver  (Rosalind: IPRB)")
    print("=" * 60)
    print(
        "\nThis program calculates the probability that two randomly\n"
        "selected organisms produce offspring with the dominant\n"
        "phenotype, given a mixed population of:\n"
        "  AA (homozygous dominant), Aa (heterozygous), aa (recessive).\n"
    )

    while True:
        population = InputHandler.get_population_from_user()

        print(f"\n{population.summary()}")

        solver = MendelSolver(population)
        result = solver.calculate()

        print(solver.detailed_breakdown())
        print(f"\nResult: {result:.5f}")

        again = input("\nSolve another? (y/n): ").strip().lower()
        if again not in ("y", "yes"):
            print("\nGoodbye.")
            break


def run_from_string(data: str, verbose: bool = True) -> float:
    """
    Solve directly from a space-separated string (e.g. Rosalind dataset format).

    Args:
        data: Space-separated string of three integers "k m n".
        verbose: If True, print the breakdown to stdout.

    Returns:
        The probability as a float.
    """
    population = InputHandler.get_population_from_string(data)
    solver = MendelSolver(population)
    result = solver.calculate()

    if verbose:
        print(population.summary())
        print(solver.detailed_breakdown())
        print(f"\nResult: {result:.5f}")

    return result


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 2:
        # Accept a quoted string argument: python mendels_first_law.py "2 2 2"
        try:
            run_from_string(sys.argv[1])
        except (ValueError, TypeError) as exc:
            print(f"[Error] {exc}")
            sys.exit(1)
    elif len(sys.argv) == 4:
        # Accept three separate arguments: python mendels_first_law.py 2 2 2
        try:
            run_from_string(" ".join(sys.argv[1:]))
        except (ValueError, TypeError) as exc:
            print(f"[Error] {exc}")
            sys.exit(1)
    else:
        # Interactive mode
        run_interactive()
