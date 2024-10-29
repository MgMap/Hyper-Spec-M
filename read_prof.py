import pstats
import argparse

def main(profile_file, top_n):
    # Load the profile data
    stats = pstats.Stats(profile_file)

    # Optionally strip directory paths for easier reading
    stats.strip_dirs()

    # Sort by cumulative time and print the top 'top_n' functions
    stats.sort_stats("cumulative").print_stats(top_n)

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Read a cProfile .prof file.')
    parser.add_argument('profile_file', type=str, help='Path to the .prof file')
    parser.add_argument('top_n', type=int, help='Number of top functions to display')

    args = parser.parse_args()
    
    main(args.profile_file, args.top_n)
