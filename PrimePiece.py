import math

class PrimePiece:
    """
    A novel tokenizer that represents sequences of prime numbers based on their
    spatial relationship on an Ulam (prime number) spiral.

    The tokenizer works by:
    1. Generating a grid representing the Ulam spiral up to a specified maximum number.
    2. Mapping each prime number to its (x, y) coordinate on this spiral.
    3. Encoding a sequence of primes as a sequence of (dx, dy) vectors,
       representing the "jump" from one prime's coordinate to the next.
    4. Decoding a sequence of vectors back into the original prime sequence,
       given a starting prime.
    """

    def __init__(self, max_n=10000):
        """
        Initializes the tokenizer and pre-computes the Ulam spiral.

        Args:
            max_n (int): The maximum number to include in the spiral.
                         A larger number allows tokenizing larger primes but
                         uses more memory and takes longer to initialize.
        """
        if max_n < 2:
            raise ValueError("max_n must be at least 2.")
        self.max_n = max_n
        self.coord_to_num = {}
        self.num_to_coord = {}
        self.primes = self._sieve(max_n)
        self.prime_coords = {}

        self._generate_spiral()
        self._map_primes_to_coords()

    def _sieve(self, n):
        """Generates a set of prime numbers up to n using a Sieve of Eratosthenes."""
        primes = [True] * (n + 1)
        if n >= 0:
            primes[0] = False
        if n >= 1:
            primes[1] = False
        for i in range(2, int(math.sqrt(n)) + 1):
            if primes[i]:
                for multiple in range(i * i, n + 1, i):
                    primes[multiple] = False
        prime_set = {i for i, is_prime in enumerate(primes) if is_prime}
        return prime_set

    def _generate_spiral(self):
        """Generates the Ulam spiral and populates the coordinate maps."""
        x, y = 0, 0
        dx, dy = 1, 0
        steps_in_leg = 1
        leg_count = 0
        turn_count = 0

        for n in range(1, self.max_n + 1):
            self.coord_to_num[(x, y)] = n
            self.num_to_coord[n] = (x, y)

            if leg_count < steps_in_leg:
                x, y = x + dx, y + dy
                leg_count += 1
            else:
                # Turn counter-clockwise
                dx, dy = -dy, dx
                leg_count = 1
                turn_count += 1
                x, y = x + dx, y + dy
                if turn_count % 2 == 0:
                    steps_in_leg += 1
        print(f"Spiral generated up to {self.max_n} with {len(self.num_to_coord)} numbers.")

    def _map_primes_to_coords(self):
        """Creates a dedicated map for prime numbers to their coordinates."""
        for prime in self.primes:
            if prime in self.num_to_coord:
                self.prime_coords[prime] = self.num_to_coord[prime]
        print(f"Found and mapped {len(self.prime_coords)} primes.")

    def encode(self, prime_sequence):
        """
        Encodes a list of prime numbers into a sequence of relational vector tokens.

        Args:
            prime_sequence (list[int]): A list of prime numbers to encode.
                                        All primes must be <= max_n.

        Returns:
            list[tuple[int, int]]: A list of (dx, dy) vector tokens.
        """
        if not prime_sequence or len(prime_sequence) < 2:
            return []

        tokens = []
        for i in range(len(prime_sequence) - 1):
            p1 = prime_sequence[i]
            p2 = prime_sequence[i+1]

            if p1 not in self.prime_coords or p2 not in self.prime_coords:
                raise ValueError(
                    f"Prime {p1} or {p2} not found in the pre-computed spiral. "
                    f"Try initializing with a larger max_n."
                )

            x1, y1 = self.prime_coords[p1]
            x2, y2 = self.prime_coords[p2]

            tokens.append((x2 - x1, y2 - y1))

        return tokens

    def decode(self, start_prime, tokens):
        """
        Decodes a sequence of vector tokens back into a list of prime numbers.

        Args:
            start_prime (int): The first prime number in the original sequence.
            tokens (list[tuple[int, int]]): The list of (dx, dy) vector tokens.

        Returns:
            list[int]: The reconstructed sequence of prime numbers.
        """
        if start_prime not in self.prime_coords:
            raise ValueError(
                f"Starting prime {start_prime} not found in the pre-computed spiral."
            )

        prime_sequence = [start_prime]
        current_x, current_y = self.prime_coords[start_prime]

        for dx, dy in tokens:
            next_x, next_y = current_x + dx, current_y + dy
            coord = (next_x, next_y)
            
            if coord not in self.coord_to_num:
                raise ValueError(f"Decoding failed: Coordinate {coord} is not on the spiral.")
            
            num = self.coord_to_num[coord]
            
            if num not in self.primes:
                print(f"Warning: Decoded number {num} at {coord} is not a prime.")

            prime_sequence.append(num)
            current_x, current_y = next_x, next_y
            
        return prime_sequence


if __name__ == '__main__':
    print("Initializing PrimePiece tokenizer...")
    # Initialize with a spiral large enough for our test primes
    tokenizer = PrimePiece(max_n=200)

    print("\n--- Encoding Example ---")
    # A sequence of the first several prime numbers
    primes_to_encode = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    print(f"Original prime sequence: {primes_to_encode}")

    # Encode the sequence into vector tokens
    encoded_tokens = tokenizer.encode(primes_to_encode)
    print(f"Encoded vector tokens:   {encoded_tokens}")
    print(f"Original sequence length: {len(primes_to_encode)}")
    print(f"Token sequence length:    {len(encoded_tokens)}")


    print("\n--- Decoding Example ---")
    # Decode the tokens back, starting from the first prime
    start_prime = primes_to_encode[0]
    decoded_primes = tokenizer.decode(start_prime, encoded_tokens)

    print(f"Tokens to decode:      {encoded_tokens}")
    print(f"Starting with prime:   {start_prime}")
    print(f"Decoded prime sequence: {decoded_primes}")

    # Verify that the decoded sequence matches the original
    assert primes_to_encode == decoded_primes
    print("\nVerification successful: Decoded sequence matches the original!")

    print("\n--- High-level concept ---")
    p1, p2 = 7, 11
    coord1 = tokenizer.prime_coords[p1]
    coord2 = tokenizer.prime_coords[p2]
    token = tokenizer.encode([p1,p2])[0]
    print(f"To get from prime {p1} at {coord1} to prime {p2} at {coord2},")
    print(f"we use the single token: {token} (dx, dy)")
