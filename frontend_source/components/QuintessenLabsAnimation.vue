<template>
  <div class="quintessen-animation-wrapper">
    <div class="animation-container" :key="animationKey">
      <svg id="morphingSVG" viewBox="0 0 1000 500" xmlns="http://www.w3.org/2000/svg">
        <!-- Real Flask SVG (Erlenmeyer) - Downloaded from Material Design Icons -->
        <!-- Morphs to Q -->
        <g transform="translate(50, 50) scale(3)">
          <path class="glassware-filled glass-1" d="M5 19a1 1 0 0 0 1 1h12a1 1 0 0 0 1-1c0-.21-.07-.41-.18-.57L13 8.35V4h-2v4.35L5.18 18.43c-.11.16-.18.36-.18.57m1 3a3 3 0 0 1-3-3c0-.6.18-1.16.5-1.63L9 7.81V6a1 1 0 0 1-1-1V4a2 2 0 0 1 2-2h4a2 2 0 0 1 2 2v1a1 1 0 0 1-1 1v1.81l5.5 9.56c.32.47.5 1.03.5 1.63a3 3 0 0 1-3 3zm7-6l1.34-1.34L16.27 18H7.73l2.66-4.61zm-.5-4a.5.5 0 0 1 .5.5a.5.5 0 0 1-.5.5a.5.5 0 0 1-.5-.5a.5.5 0 0 1 .5-.5">
            <animate
              attributeName="d"
              dur="1.5s"
              begin="2s"
              fill="freeze"
              values="
                M5 19a1 1 0 0 0 1 1h12a1 1 0 0 0 1-1c0-.21-.07-.41-.18-.57L13 8.35V4h-2v4.35L5.18 18.43c-.11.16-.18.36-.18.57m1 3a3 3 0 0 1-3-3c0-.6.18-1.16.5-1.63L9 7.81V6a1 1 0 0 1-1-1V4a2 2 0 0 1 2-2h4a2 2 0 0 1 2 2v1a1 1 0 0 1-1 1v1.81l5.5 9.56c.32.47.5 1.03.5 1.63a3 3 0 0 1-3 3zm7-6l1.34-1.34L16.27 18H7.73l2.66-4.61zm-.5-4a.5.5 0 0 1 .5.5a.5.5 0 0 1-.5.5a.5.5 0 0 1-.5-.5a.5.5 0 0 1 .5-.5;
                M 12 6 Q 6 6 6 12 Q 6 18 12 18 Q 18 18 18 12 Q 18 6 12 6 M 16 16 L 20 20;
                M 12 6 Q 6 6 6 12 Q 6 18 12 18 Q 18 18 18 12 Q 18 6 12 6 M 16 16 L 20 20
              "
            />
            <animate attributeName="opacity" dur="0.5s" begin="2s" fill="freeze" values="1; 1; 0" />
          </path>
        </g>

        <!-- Real Test Tube SVG - Downloaded from Material Design Icons -->
        <g transform="translate(250, 280) scale(3)">
          <path class="glassware-filled glass-2" d="M7 2v2h1v14a4 4 0 0 0 4 4a4 4 0 0 0 4-4V4h1V2zm4 14c-.6 0-1-.4-1-1s.4-1 1-1s1 .4 1 1s-.4 1-1 1m2-4c-.6 0-1-.4-1-1s.4-1 1-1s1 .4 1 1s-.4 1-1 1m1-5h-4V4h4z">
            <animate
              attributeName="d"
              dur="1.5s"
              begin="2.2s"
              fill="freeze"
              values="
                M7 2v2h1v14a4 4 0 0 0 4 4a4 4 0 0 0 4-4V4h1V2zm4 14c-.6 0-1-.4-1-1s.4-1 1-1s1 .4 1 1s-.4 1-1 1m2-4c-.6 0-1-.4-1-1s.4-1 1-1s1 .4 1 1s-.4 1-1 1m1-5h-4V4h4z;
                M7 4v16h2V4z M7 20h10;
                M0 0v20h4v0z M0 20h16
              "
            />
            <animate attributeName="opacity" dur="0.5s" begin="2.2s" fill="freeze" values="1; 1; 0" />
          </path>
        </g>

        <!-- Real Beaker SVG - Downloaded from Material Design Icons -->
        <g transform="translate(420, 280) scale(3)">
          <path class="glassware-filled glass-3" d="M3 3h18v2a2 2 0 0 0-2 2v12a2 2 0 0 1-2 2H7a2 2 0 0 1-2-2V7a2 2 0 0 0-2-2zm4 2v2h5v1H7v1h3v1H7v1h3v1H7v1h5v1H7v1h3v1H7v3h10V5z">
            <animate
              attributeName="d"
              dur="1.5s"
              begin="2.4s"
              fill="freeze"
              values="
                M3 3h18v2a2 2 0 0 0-2 2v12a2 2 0 0 1-2 2H7a2 2 0 0 1-2-2V7a2 2 0 0 0-2-2zm4 2v2h5v1H7v1h3v1H7v1h3v1H7v1h5v1H7v1h3v1H7v3h10V5z;
                M4 4h14v14H4z M6 8h10;
                M0 20L12 0L24 20 M6 12h12
              "
            />
            <animate attributeName="opacity" dur="0.5s" begin="2.4s" fill="freeze" values="1; 1; 0" />
          </path>
        </g>

        <!-- Real Bottle SVG - Downloaded from Material Design Icons -->
        <g transform="translate(570, 280) scale(3)">
          <path class="glassware-filled glass-4" d="M13 4h-2l-1-2h4zm6 9v9H5v-9c0-2.76 2.24-5 5-5V6H9V5h6v1h-1v2c2.76 0 5 2.24 5 5">
            <animate
              attributeName="d"
              dur="1.5s"
              begin="2.6s"
              fill="freeze"
              values="
                M13 4h-2l-1-2h4zm6 9v9H5v-9c0-2.76 2.24-5 5-5V6H9V5h6v1h-1v2c2.76 0 5 2.24 5 5;
                M10 2h4v20H10z M10 2h8 M10 11h8;
                M0 0v20h4 M0 0h8 Q12 0 12 4 Q12 8 8 8h-4 M4 8h4 Q12 8 12 12 Q12 20 8 20 h-8
              "
            />
            <animate attributeName="opacity" dur="0.5s" begin="2.6s" fill="freeze" values="1; 1; 0" />
          </path>
        </g>

        <!-- Another test tube for S -->
        <g transform="translate(720, 280) scale(3)">
          <path class="glassware-filled glass-4" style="animation-delay: 1s;" d="M7 2v2h1v14a4 4 0 0 0 4 4a4 4 0 0 0 4-4V4h1V2z" opacity="0">
            <animate
              attributeName="d"
              dur="1.2s"
              begin="2.8s"
              fill="freeze"
              values="
                M7 2v2h1v14a4 4 0 0 0 4 4a4 4 0 0 0 4-4V4h1V2z;
                M8 0 Q4 0 4 4 Q4 7 8 7 Q12 7 12 10 Q12 14 8 14 Q4 14 4 16 Q4 20 8 20
              "
            />
            <animate attributeName="opacity" dur="0.3s" begin="2.8s" fill="freeze" values="0; 1; 0" />
          </path>
        </g>

        <!-- Sans-serif Q (final form) -->
        <text x="110" y="210" font-family="Arial, Helvetica, sans-serif" font-size="120" font-weight="bold" class="text-el">
          Q
          <animate attributeName="opacity" dur="0.6s" begin="3.5s" fill="freeze" values="0; 1" />
        </text>

        <!-- "uintessen" slides out from Q -->
        <text x="240" y="210" font-family="Arial, Helvetica, sans-serif" font-size="90" font-weight="300" class="text-el" opacity="0">
          uintessen
          <animate attributeName="opacity" dur="0.8s" begin="4s" fill="freeze" values="0; 1" />
          <animateTransform
            attributeName="transform"
            type="translate"
            dur="0.8s"
            begin="4s"
            fill="freeze"
            values="-80,0; 0,0"
          />
        </text>

        <!-- LABS (final form) -->
        <text x="110" y="390" font-family="Arial, Helvetica, sans-serif" font-size="90" font-weight="300" letter-spacing="15" class="text-el">
          LABS
          <animate attributeName="opacity" dur="0.6s" begin="3.5s" fill="freeze" values="0; 1" />
        </text>
      </svg>
    </div>

    <q-btn
      v-if="showRestartButton"
      class="restart-btn"
      label="Restart Animation"
      @click="restartAnimation"
      flat
      color="dark"
    />
  </div>
</template>

<script setup>
import { ref, onMounted, onUnmounted } from 'vue'

// Props
const props = defineProps({
  autoRestart: {
    type: Boolean,
    default: true
  },
  restartInterval: {
    type: Number,
    default: 9000
  },
  showRestartButton: {
    type: Boolean,
    default: true
  }
})

// State
const animationKey = ref(0)
let intervalId = null

// Methods
const restartAnimation = () => {
  animationKey.value++
}

// Lifecycle hooks
onMounted(() => {
  if (props.autoRestart) {
    intervalId = setInterval(() => {
      restartAnimation()
    }, props.restartInterval)
  }
})

onUnmounted(() => {
  if (intervalId) {
    clearInterval(intervalId)
  }
})

// Expose methods for parent component
defineExpose({
  restartAnimation
})
</script>

<style scoped>
.quintessen-animation-wrapper {
  width: 100%;
  min-height: 500px;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  position: relative;
  background: #ffffff;
}

.animation-container {
  position: relative;
  width: 1000px;
  height: 500px;
  max-width: 100%;
}

#morphingSVG {
  width: 100%;
  height: 100%;
}

/* Real lab glassware icons */
.glassware-filled {
  fill: #000000;
  opacity: 0;
  animation: fadeInGlass 1s ease-out forwards;
}

.glass-1 { animation-delay: 0.2s; }
.glass-2 { animation-delay: 0.4s; }
.glass-3 { animation-delay: 0.6s; }
.glass-4 { animation-delay: 0.8s; }

/* Sans-serif text */
.text-el {
  fill: #000000;
  opacity: 0;
}

@keyframes fadeInGlass {
  to {
    opacity: 1;
  }
}

/* Restart button */
.restart-btn {
  margin-top: 40px;
  border: 1px solid #000000;
}

/* Responsive design */
@media (max-width: 768px) {
  .animation-container {
    width: 100%;
    max-width: 700px;
    height: 400px;
  }

  #morphingSVG {
    transform: scale(0.8);
  }
}

@media (max-width: 480px) {
  .animation-container {
    width: 100%;
    height: 350px;
  }

  #morphingSVG {
    transform: scale(0.6);
  }
}
</style>
