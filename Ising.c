#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef int8_t spin_t;
const spin_t SPIN_UP = 1;
const spin_t SPIN_DOWN = -1;

typedef int64_t energy_t;
typedef double weight_t;

weight_t uniformRandom() {
    return ((weight_t) rand() / (weight_t) RAND_MAX);
}

typedef struct {
    size_t width;
    size_t height;
    spin_t* data;
} Lattice;

Lattice allocateLattice(size_t width, size_t height) {
    Lattice result;
    result.width = width;
    result.height = height;
    result.data = calloc(width * height, sizeof(spin_t));
    return result;
}

void deallocateLattice(Lattice l) {
    free(l.data);
}

size_t getIndex(Lattice lattice, size_t x, size_t y) {
    x = (x + lattice.width) % lattice.width;
    y = (y + lattice.height) % lattice.height;
    return x * lattice.height + y;
}

spin_t get(Lattice lattice, size_t x, size_t y) {
    return lattice.data[getIndex(lattice, x, y)];
}

spin_t set(Lattice lattice, size_t x, size_t y, spin_t value) {
    lattice.data[getIndex(lattice, x, y)] = value;
}

void fillGroundState(Lattice lattice, spin_t spin) {
    size_t x, y;
    for (x = 0; x < lattice.width; ++x)
        for (y = 0; y < lattice.height; ++y)
            set(lattice, x, y, spin);
}

void flip(Lattice lattice, size_t x, size_t y) {
    set(lattice, x, y, - get(lattice, x, y));
}

energy_t calcEnergy(Lattice lattice) {
    size_t x, y;
    energy_t energy = 0;
    for (x = 0; x < lattice.width; ++x) {
        for (y = 0; y < lattice.height; ++y) {
            energy += get(lattice, x, y) * get(lattice, x + 1, y);
            energy += get(lattice, x, y) * get(lattice, x - 1, y);
            energy += get(lattice, x, y) * get(lattice, x, y + 1);
            energy += get(lattice, x, y) * get(lattice, x, y - 1);
        }
    }
    return energy;
}

weight_t calcWeight(Lattice lattice, weight_t beta) {
    energy_t energy = calcEnergy(lattice);
    return exp(- beta * energy);
}

void tryFlipOneSpin(Lattice lattice, weight_t beta, weight_t* buffer) {
    size_t x, y;
    for (x = 0; x < lattice.width; ++x) {
        for (y = 0; y < lattice.height; ++y) {
            flip(lattice, x, y);
            weight_t weight = calcWeight(lattice, beta);
            buffer[getIndex(lattice, x, y)] = weight;
            flip(lattice, x, y);
        }
    }

    weight_t sum = 0;
    for (x = 0; x < lattice.width; ++x)
        for (y = 0; y < lattice.height; ++y)
            sum += buffer[getIndex(lattice, x, y)];

    for (x = 0; x < lattice.width; ++x)
        for (y = 0; y < lattice.height; ++y)
            buffer[getIndex(lattice, x, y)] /= sum;

    double r = uniformRandom();

    sum = 0;
    for (x = 0; x < lattice.width; ++x) {
        for (y = 0; y < lattice.height; ++y) {
            sum += buffer[getIndex(lattice, x, y)];
            if (sum > r) {
                flip(lattice, x, y);
                return;
            }
        }
    }
}

void evolveIntoThermalState(Lattice lattice, weight_t beta, size_t max_consequential_hits) {
    energy_t energy = calcEnergy(lattice);
    weight_t *buffer = calloc(lattice.width * lattice.height, sizeof(weight_t));

    size_t consequential_hits = 0;

    for (;;) {
        tryFlipOneSpin(lattice, beta, buffer);
        energy_t newEnergy = calcEnergy(lattice);

        if (energy == newEnergy)
            ++consequential_hits;
        else
            consequential_hits = 0;

        if (consequential_hits >= max_consequential_hits)
            break;

        energy = newEnergy;
    }

    free(buffer);
}

void printLattice(Lattice lattice, FILE* file_descriptor) {
    size_t x, y;
    for (x = 0; x < lattice.width; ++x) {
        for (y = 0; y < lattice.height; ++y)
            fprintf(file_descriptor, "%c", get(lattice, x, y) == 1 ? '+' : '-');
        fprintf(file_descriptor, "\n");
    }
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    Lattice lattice = allocateLattice(5, 5);

    fillGroundState(lattice, SPIN_DOWN);
    printf("Lattice at the ground state:\n");
    printLattice(lattice, stdout);

    energy_t groundEnergy = calcEnergy(lattice);
    printf("Energy: %lld\n", groundEnergy);

    evolveIntoThermalState(lattice, 1, 3);
    printf("Lattice at the thermal state:\n");
    printLattice(lattice, stdout);

    energy_t thermalEnergy = calcEnergy(lattice);
    printf("Energy: %lld\n", thermalEnergy);

    deallocateLattice(lattice);

    return 0;
}
