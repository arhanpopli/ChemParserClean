# QuintessenLabs Animation

A stunning chemistry-themed animation featuring beakers transforming into the QuintessenLabs logo.

## Features

- **Chemistry-themed**: Beautiful SVG beakers with animated liquid
- **Smooth transformations**: Beakers morph into the logo text
- **Purple branding**: "Quintessen" appears in vibrant purple (#a855f7)
- **Auto-restart**: Optional automatic animation looping
- **Responsive**: Works on desktop and mobile devices
- **Multiple formats**: Standalone HTML and Vue component versions

## Animation Sequence

1. **0.2s - 0.8s**: Four chemistry beakers fade in sequentially from left to right
2. **1.5s - 2.5s**: Beakers transform - first beaker morphs into "Q", others into "LABS"
3. **2.5s - 3s**: The "Q" letter appears with rotation effect, "LABS" materializes below
4. **3s - 4s**: "uintessen" slides out from the Q, forming complete "Quintessen" text

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
- **CSS3 Animations**: @keyframes, transforms, transitions
- **SVG**: Scalable vector graphics for beakers

### Browser Compatibility

- ✅ Chrome/Edge (latest)
- ✅ Firefox (latest)
- ✅ Safari (latest)
- ✅ Mobile browsers (iOS Safari, Chrome Mobile)

### Performance

- GPU-accelerated CSS animations
- No external animation libraries required
- Lightweight SVG graphics
- Optimized for 60fps on modern devices

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
