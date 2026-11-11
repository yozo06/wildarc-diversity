# @wildarc/diversity

Biodiversity indices from species abundances: richness, Shannon, Simpson, Pielou's evenness and Hill numbers. It is
one pure function with no dependencies, and every formula is cited. When a measure does not apply to the input, the
result says so with a reason; it never returns `NaN`, `Infinity` or an invented fallback number.

## Install

The package is not on npm yet. Install it from this repository, pinned to a tag:

```bash
npm install github:yozo06/wildarc-diversity#v0.1.0
```

npm compiles the package while installing it (the `prepare` script builds `dist/`), so the first install is slower
than one from the registry. What arrives is CommonJS with type declarations and no runtime dependencies.

## Usage

```ts
import { diversityProfile } from '@wildarc/diversity';

const profile = diversityProfile({
    'Grevillea robusta': 50,
    'Coffea arabica': 30,
    'Mangifera indica': 20,
});
```

```json
{
  "totalAbundance": 100,
  "richness": 3,
  "shannon": { "value": 1.0296530140645737 },
  "pielouEvenness": { "value": 0.9372305632161297 },
  "simpsonLambda": { "value": 0.3737373737373737 },
  "giniSimpson": { "value": 0.6262626262626263 },
  "hill": {
    "q0": 3,
    "q1": { "value": 2.800094072853832 },
    "q2": { "value": 2.6315789473684212 }
  }
}
```

The input is either a record keyed by species name or a plain array of abundances (`diversityProfile([50, 30, 20])`
gives the same result). A species with abundance 0 is absent: it counts towards neither richness nor the total.

## Measures

`nᵢ` is the abundance of species *i*, `N` the total, `S` the number of species present and `pᵢ = nᵢ / N`.

| Field            | Formula                    | Source                                  |
| ---------------- | -------------------------- | --------------------------------------- |
| `totalAbundance` | N = Σ nᵢ                   |                                         |
| `richness`       | S                          | Magurran (2004), ch. 4                  |
| `shannon`        | H′ = −Σ pᵢ ln pᵢ           | Magurran (2004), ch. 4                  |
| `pielouEvenness` | J = H′ / ln S              | Magurran (2004), ch. 4                  |
| `simpsonLambda`  | λ = Σ nᵢ(nᵢ−1) / (N(N−1))  | Magurran (2004), ch. 4                  |
| `giniSimpson`    | 1 − λ                      | Hurlbert (1971), as PIE; see Notes      |
| `hill.q0`        | S                          | Jost (2006)                             |
| `hill.q1`        | exp(H′)                    | Jost (2006)                             |
| `hill.q2`        | 1 / Σ pᵢ²                  | Jost (2006)                             |

- Magurran, A. E. (2004). *Measuring Biological Diversity*. Blackwell Publishing.
- Jost, L. (2006). Entropy and diversity. *Oikos* 113(2):363–375.
- Hurlbert, S. H. (1971). The nonconcept of species diversity: a critique and alternative parameters. *Ecology*
  52(4):577–586.

Every measure except `totalAbundance`, `richness` and `hill.q0` is a `Measure`: either `{ value: number }` or
`{ value: null, reason: string }`.

## When a measure is refused

| `reason`                                             | Fields                                   | When                                 |
| ---------------------------------------------------- | ---------------------------------------- | ------------------------------------ |
| `empty community: no individuals counted`            | all six measures                         | no abundance is greater than 0       |
| `evenness is undefined for fewer than 2 species`     | `pielouEvenness`                         | S = 1 (ln 1 = 0, so J would be 0/0)  |
| `Simpson's unbiased index needs whole-number counts` | `simpsonLambda`, `giniSimpson`           | any abundance is not a whole number  |
| `Simpson's index needs at least 2 individuals`       | `simpsonLambda`, `giniSimpson`           | N < 2                                |

The whole-number check comes first, so `[0.5]` reports the whole-number reason rather than the N < 2 one.

## Invalid input

`diversityProfile` throws a `RangeError` when:

- an abundance is negative, `NaN` or infinite (the message names the species or array index);
- the abundances sum to more than the largest finite number;
- `logBase` is not a finite number greater than 0, or is exactly 1.

## Notes

**`hill.q2` is not `1 / simpsonLambda`.** `hill.q2` is the classic inverse Simpson index, calculated from the
proportions (the "plug-in" form, as if drawing with replacement). `simpsonLambda` is the unbiased finite-sample form,
which draws without replacement. In the example above `1 / simpsonLambda` is about 2.676 while `hill.q2` is about
2.632. The two converge as N grows.

**`giniSimpson` is the unbiased form, not the one in Jost (2006).** `giniSimpson` is 1 − λ with the unbiased λ: the
chance that two individuals drawn without replacement belong to different species. Written out, it is
(N/(N−1))(1 − Σ pᵢ²), which Hurlbert (1971) called the probability of interspecific encounter (PIE). The
Gini-Simpson index in Jost (2006), Table 1 is the plug-in form 1 − Σ pᵢ². To get that one, use `1 - 1 / hill.q2`.
In the example above `giniSimpson` is about 0.626 and the plug-in form is 0.62. The gap is widest for two
individuals of two species: 1 against 0.5.

**Log bases.** H′ in base *b* is H′ in natural logs divided by ln *b*. The `logBase` option does this for you:
`diversityProfile([10, 10, 10, 10], { logBase: 2 }).shannon` is 2 bits, where the default gives ln 4 ≈ 1.386 nats.
`logBase` affects `shannon` only. Evenness and the Hill numbers do not depend on the base.

**Fractional abundances** (cover, biomass, basal area) are accepted. Shannon, evenness and the Hill numbers are
calculated from the proportions; the unbiased Simpson index is refused because it is defined for counts of individuals.

**Floating point.** Results are IEEE 754 doubles and are not clamped. A perfectly even community can report an
evenness or Hill number one unit in the last place past its mathematical bound (for 5 equally abundant species,
J = 1.0000000000000002). Compare with a tolerance rather than `<= 1`.

**Very small shares.** A species whose share of the total is too small to represent (below about 1e-308) still counts
towards richness but adds nothing to Shannon or the Hill numbers.

## Develop

```bash
npm install        # also builds dist/ through the prepare script
npm test           # 28 tests, including worked examples from the sources
npm run typecheck
npm run build
```
