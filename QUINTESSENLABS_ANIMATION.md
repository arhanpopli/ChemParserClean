# QuintessenLabs Animation

A stunning chemistry-themed animation featuring beakers transforming into the QuintessenLabs logo.

## Features

- **Line morphing animation**: Beaker outlines literally transform and bend into text
- **Different sized beakers**: 4 beakers of varying sizes (small, medium, large, medium-small)
- **SVG path morphing**: Lines smoothly morph from beaker shapes into letter shapes
- **Purple branding**: "Quintessen" appears in vibrant purple (#a855f7)
- **Line drawing effect**: Uses stroke-dasharray for smooth line drawing animations
- **Auto-restart**: Optional automatic animation looping
- **Responsive**: Works on desktop and mobile devices
- **Multiple formats**: Standalone HTML and Vue component versions

## Animation Sequence

1. **0.2s - 1.7s**: Four beakers of different sizes draw in sequentially (line drawing effect)
2. **2s - 3.5s**: Beaker outlines morph and bend into letters:
   - Small beaker → "Q" (with tail)
   - Medium beaker → "u"
   - Large beaker → "i" (with dot)
   - Medium-small beaker → "n"
3. **2.8s - 3.4s**: "LABS" letters draw in below (L, A, B, S)
4. **3.5s - 5s**: Remaining letters of "Quintessen" draw in (t, e, s, s, e, n)

## Files Created

```
ChemParserClean/
├── quintessenlabs-animation.html              # Standalone HTML version
├── frontend_source/
│   ├── components/
│   │   └── QuintessenLabsAnimation.vue       # Vue component
│   └── pages/
│       └── QuintessenLabsDemo.vue            # Demo page with controls
└── QUINTESSENLABS_ANIMATION.md               # This documentation
```

## Usage

### Option 1: Standalone HTML

Simply open `quintessenlabs-animation.html` in any web browser:

```bash
# Open directly in browser
open quintessenlabs-animation.html

# Or serve with a local server
python -m http.server 8000
# Then visit: http://localhost:8000/quintessenlabs-animation.html
```

### Option 2: Vue Component (Recommended for App Integration)

#### Basic Usage

```vue
<template>
  <QuintessenLabsAnimation />
</template>

<script setup>
import QuintessenLabsAnimation from 'components/QuintessenLabsAnimation.vue'
</script>
```

#### Advanced Usage with Props

```vue
<template>
  <QuintessenLabsAnimation
    ref="animationRef"
    :auto-restart="true"
    :restart-interval="10000"
    :show-restart-button="false"
  />

  <q-btn @click="triggerAnimation">
    Restart Animation
  </q-btn>
</template>

<script setup>
import { ref } from 'vue'
import QuintessenLabsAnimation from 'components/QuintessenLabsAnimation.vue'

const animationRef = ref(null)

const triggerAnimation = () => {
  animationRef.value?.restartAnimation()
}
</script>
```

#### Props

| Prop | Type | Default | Description |
|------|------|---------|-------------|
| `autoRestart` | Boolean | `true` | Automatically restart the animation |
| `restartInterval` | Number | `8000` | Time in milliseconds between auto-restarts |
| `showRestartButton` | Boolean | `true` | Show the restart button below animation |

#### Exposed Methods

- `restartAnimation()`: Manually trigger the animation to restart

### Option 3: Demo Page

To see the animation with interactive controls, add the demo page to your router:

1. **Update `frontend_source/router/routes.js`**:

```javascript
const routes = [
  // ... existing routes
  {
    path: '/quintessen-demo',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      {
        path: '',
        component: () => import('pages/QuintessenLabsDemo.vue')
      }
    ]
  }
]
```

2. **Visit the demo page**: Navigate to `/quintessen-demo` in your app

## Customization

### Colors

The animation uses the following color scheme:

- **Purple (Quintessen)**: `#a855f7` - Main brand color
- **Gray (Labs)**: `#64748b` - Secondary text color
- **Beaker gradient**: `#e0e7ff` to `#c7d2fe` - Light purple gradient
- **Background**: `linear-gradient(135deg, #1a1a2e 0%, #16213e 100%)`

To customize colors, modify the CSS variables in the component or HTML file:

```css
.letter-q, .quintessen-text {
  color: #a855f7; /* Change to your brand color */
}

.labs-text {
  color: #64748b; /* Change to your secondary color */
}
```

### Animation Timing

Adjust animation speed by modifying the animation durations and delays:

```css
/* Fade in duration */
@keyframes fadeIn {
  /* ... 0.6s can be changed */
}

/* Delay between beakers appearing */
.beaker:nth-child(1) {
  animation-delay: 0.2s; /* Adjust delays */
}
```

### Size

For the Vue component, adjust the container size:

```css
.animation-container {
  width: 600px;  /* Adjust width */
  height: 300px; /* Adjust height */
}
```

## Technical Details

### Technologies Used

- **Vue 3**: Composition API with `<script setup>`
- **Quasar Framework**: For buttons and UI components (Vue version only)
- **SVG SMIL Animations**: Native SVG `<animate>` elements for path morphing
- **CSS3**: stroke-dasharray/stroke-dashoffset for line drawing effects
- **SVG Path Commands**: Complex path morphing with L (line), Q (quadratic curve), M (move) commands

### Browser Compatibility

- ✅ Chrome/Edge (latest) - Full support for SMIL animations
- ✅ Firefox (latest) - Full support for SMIL animations
- ✅ Safari (latest) - Full support for SMIL animations
- ✅ Mobile browsers (iOS Safari, Chrome Mobile)

### Performance

- GPU-accelerated SVG animations
- No external animation libraries required (uses native SVG SMIL)
- Lightweight stroke-only graphics (no fills)
- Optimized for 60fps on modern devices
- Smooth path morphing with cubic bezier timing functions

## Examples

### Landing Page

```vue
<template>
  <div class="landing-page">
    <QuintessenLabsAnimation
      :auto-restart="true"
      :restart-interval="12000"
      :show-restart-button="false"
    />
    <h1>Welcome to QuintessenLabs</h1>
  </div>
</template>
```

### Loading Screen

```vue
<template>
  <div class="loading-screen">
    <QuintessenLabsAnimation
      :auto-restart="false"
      :show-restart-button="false"
    />
    <p>Loading your chemistry tools...</p>
  </div>
</template>
```

### Modal/Dialog

```vue
<template>
  <q-dialog v-model="showWelcome">
    <q-card>
      <q-card-section>
        <QuintessenLabsAnimation
          :auto-restart="false"
          :show-restart-button="false"
        />
      </q-card-section>
      <q-card-section>
        <h2>Welcome to QuintessenLabs!</h2>
        <p>Your chemistry companion</p>
      </q-card-section>
    </q-card>
  </q-dialog>
</template>
```

## Troubleshooting

### Animation doesn't restart automatically

Make sure `autoRestart` prop is set to `true`:

```vue
<QuintessenLabsAnimation :auto-restart="true" />
```

### Animation appears cut off

Ensure the parent container has enough height:

```css
.parent-container {
  min-height: 300px;
}
```

### Colors don't match brand

Update the color values in the `<style>` section of the component to match your brand colors.

## Credits

Created for the ChemParserClean project - A chemistry to LaTeX conversion system.

## License

Part of the ChemParserClean project.
