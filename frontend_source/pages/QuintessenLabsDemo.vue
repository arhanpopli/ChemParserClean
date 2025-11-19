<template>
  <q-page class="demo-page">
    <div class="page-container">
      <div class="header-section">
        <h1 class="page-title">QuintessenLabs Animation Demo</h1>
        <p class="page-description">
          A stunning chemistry-themed animation featuring beakers transforming into the QuintessenLabs logo
        </p>
      </div>

      <!-- Main Animation Display -->
      <q-card class="animation-card">
        <q-card-section class="animation-section">
          <QuintessenLabsAnimation
            ref="animationRef"
            :auto-restart="autoRestart"
            :restart-interval="restartInterval"
            :show-restart-button="showRestartButton"
          />
        </q-card-section>
      </q-card>

      <!-- Controls -->
      <q-card class="controls-card">
        <q-card-section>
          <h2 class="controls-title">Animation Controls</h2>

          <div class="controls-grid">
            <q-toggle
              v-model="autoRestart"
              label="Auto-restart animation"
              color="primary"
            />

            <q-toggle
              v-model="showRestartButton"
              label="Show restart button"
              color="primary"
            />

            <div class="control-item">
              <label class="control-label">Restart Interval (seconds)</label>
              <q-slider
                v-model="restartIntervalSeconds"
                :min="3"
                :max="15"
                :step="1"
                label
                color="primary"
              />
            </div>

            <q-btn
              label="Trigger Animation Now"
              color="primary"
              @click="triggerAnimation"
              class="trigger-btn"
            />
          </div>
        </q-card-section>
      </q-card>

      <!-- Code Examples -->
      <q-card class="code-card">
        <q-card-section>
          <h2 class="code-title">Usage Example</h2>

          <q-tabs v-model="selectedTab" align="left" class="code-tabs">
            <q-tab name="basic" label="Basic Usage" />
            <q-tab name="advanced" label="Advanced" />
            <q-tab name="props" label="Props" />
          </q-tabs>

          <q-tab-panels v-model="selectedTab" class="code-panels">
            <q-tab-panel name="basic">
              <pre class="code-block"><code>&lt;template&gt;
  &lt;QuintessenLabsAnimation /&gt;
&lt;/template&gt;

&lt;script setup&gt;
import QuintessenLabsAnimation from 'components/QuintessenLabsAnimation.vue'
&lt;/script&gt;</code></pre>
            </q-tab-panel>

            <q-tab-panel name="advanced">
              <pre class="code-block"><code>&lt;template&gt;
  &lt;QuintessenLabsAnimation
    ref="animationRef"
    :auto-restart="true"
    :restart-interval="10000"
    :show-restart-button="false"
  /&gt;

  &lt;q-btn @click="triggerAnimation"&gt;
    Restart Animation
  &lt;/q-btn&gt;
&lt;/template&gt;

&lt;script setup&gt;
import { ref } from 'vue'
import QuintessenLabsAnimation from 'components/QuintessenLabsAnimation.vue'

const animationRef = ref(null)

const triggerAnimation = () =&gt; {
  animationRef.value?.restartAnimation()
}
&lt;/script&gt;</code></pre>
            </q-tab-panel>

            <q-tab-panel name="props">
              <div class="props-table">
                <table>
                  <thead>
                    <tr>
                      <th>Prop</th>
                      <th>Type</th>
                      <th>Default</th>
                      <th>Description</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      <td><code>autoRestart</code></td>
                      <td>Boolean</td>
                      <td>true</td>
                      <td>Automatically restart the animation</td>
                    </tr>
                    <tr>
                      <td><code>restartInterval</code></td>
                      <td>Number</td>
                      <td>8000</td>
                      <td>Time in milliseconds between auto-restarts</td>
                    </tr>
                    <tr>
                      <td><code>showRestartButton</code></td>
                      <td>Boolean</td>
                      <td>true</td>
                      <td>Show the restart button below animation</td>
                    </tr>
                  </tbody>
                </table>
              </div>
            </q-tab-panel>
          </q-tab-panels>
        </q-card-section>
      </q-card>

      <!-- Animation Details -->
      <q-card class="details-card">
        <q-card-section>
          <h2 class="details-title">Animation Sequence</h2>
          <q-timeline color="primary">
            <q-timeline-entry
              title="Step 1: Beakers Appear"
              subtitle="0.2s - 0.8s"
              icon="science"
            >
              Four chemistry beakers fade in sequentially from left to right
            </q-timeline-entry>

            <q-timeline-entry
              title="Step 2: Transformation Begins"
              subtitle="1.5s - 2.5s"
              icon="auto_fix_high"
            >
              The first beaker morphs into the letter "Q" while the other three beakers transform into "LABS"
            </q-timeline-entry>

            <q-timeline-entry
              title="Step 3: Text Reveal"
              subtitle="2.5s - 3s"
              icon="text_fields"
            >
              The "Q" appears with a rotation effect, and "LABS" materializes below
            </q-timeline-entry>

            <q-timeline-entry
              title="Step 4: Complete Logo"
              subtitle="3s - 4s"
              icon="check_circle"
            >
              "uintessen" slides out from the Q, forming the complete "Quintessen" text in purple
            </q-timeline-entry>
          </q-timeline>
        </q-card-section>
      </q-card>
    </div>
  </q-page>
</template>

<script setup>
import { ref, computed } from 'vue'
import QuintessenLabsAnimation from 'components/QuintessenLabsAnimation.vue'

// State
const animationRef = ref(null)
const autoRestart = ref(true)
const showRestartButton = ref(true)
const restartIntervalSeconds = ref(8)
const selectedTab = ref('basic')

// Computed
const restartInterval = computed(() => restartIntervalSeconds.value * 1000)

// Methods
const triggerAnimation = () => {
  animationRef.value?.restartAnimation()
}
</script>

<style scoped>
.demo-page {
  background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
  min-height: 100vh;
  padding: 40px 20px;
}

.page-container {
  max-width: 1200px;
  margin: 0 auto;
}

.header-section {
  text-align: center;
  margin-bottom: 40px;
  color: white;
}

.page-title {
  font-size: 48px;
  font-weight: 700;
  margin-bottom: 16px;
  background: linear-gradient(135deg, #a855f7 0%, #667eea 100%);
  -webkit-background-clip: text;
  -webkit-text-fill-color: transparent;
  background-clip: text;
}

.page-description {
  font-size: 18px;
  color: #94a3b8;
  max-width: 600px;
  margin: 0 auto;
}

.animation-card {
  background: rgba(255, 255, 255, 0.05);
  backdrop-filter: blur(10px);
  border: 1px solid rgba(255, 255, 255, 0.1);
  margin-bottom: 30px;
}

.animation-section {
  padding: 60px 20px;
  background: rgba(0, 0, 0, 0.2);
}

.controls-card {
  background: rgba(255, 255, 255, 0.95);
  margin-bottom: 30px;
}

.controls-title {
  font-size: 24px;
  font-weight: 600;
  margin-bottom: 24px;
  color: #1e293b;
}

.controls-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
  gap: 24px;
  align-items: center;
}

.control-item {
  display: flex;
  flex-direction: column;
  gap: 8px;
}

.control-label {
  font-size: 14px;
  font-weight: 500;
  color: #475569;
}

.trigger-btn {
  width: 100%;
}

.code-card {
  background: rgba(255, 255, 255, 0.95);
  margin-bottom: 30px;
}

.code-title {
  font-size: 24px;
  font-weight: 600;
  margin-bottom: 24px;
  color: #1e293b;
}

.code-tabs {
  margin-bottom: 16px;
}

.code-block {
  background: #1e293b;
  color: #e2e8f0;
  padding: 20px;
  border-radius: 8px;
  overflow-x: auto;
  font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
  font-size: 14px;
  line-height: 1.6;
}

.code-block code {
  color: inherit;
}

.props-table {
  overflow-x: auto;
}

.props-table table {
  width: 100%;
  border-collapse: collapse;
}

.props-table th,
.props-table td {
  padding: 12px;
  text-align: left;
  border-bottom: 1px solid #e2e8f0;
}

.props-table th {
  background: #f8fafc;
  font-weight: 600;
  color: #1e293b;
}

.props-table code {
  background: #f1f5f9;
  padding: 2px 6px;
  border-radius: 4px;
  font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
  font-size: 13px;
  color: #a855f7;
}

.details-card {
  background: rgba(255, 255, 255, 0.95);
  margin-bottom: 30px;
}

.details-title {
  font-size: 24px;
  font-weight: 600;
  margin-bottom: 24px;
  color: #1e293b;
}

/* Responsive */
@media (max-width: 768px) {
  .page-title {
    font-size: 32px;
  }

  .page-description {
    font-size: 16px;
  }

  .controls-grid {
    grid-template-columns: 1fr;
  }
}
</style>
